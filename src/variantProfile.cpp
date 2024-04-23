#include "variantProfile.hpp"
#include "eigen3/Eigen/src/SparseCore/SparseMatrix.h"
#include "genotypeDistribution.hpp"

VariantProfile::VariantProfile()
{
    this->sMinMapped = 1;
    this->sMaxMapped = 1000;
    this->readLength = 150;
    this->overlap = 20;
    this->variantPresent = false;
    this->name = "unnamed";
}

VariantProfile::VariantProfile(std::string filename)
{
    readProfile(filename);
}

VariantProfile::VariantProfile(
    complexVariant variant,
    int filterMargin,
    int overlap,
    int readLength, 
    int sMin,
    int sMax,
    const ContigInfo & contigInfos
) : variant(variant), filterMargin(filterMargin), overlap(overlap), readLength(readLength), sMin(sMin), sMax(sMax), cInfo(contigInfos)
{    
    this->variant.createAlleleMaps(this->filterMargin, cInfo.getContigLengths());
    this->filter = ReadPairFilter(this->variant.getAllBreakpoints(), this->filterMargin, 0);
    this->sMinMapped = sMin;
    this->sMaxMapped = sMax;
    this->variantPresent = true;
    this->name = this->variant.getName();
    determinePossibleGroups();
    initMasks();
}

void VariantProfile::determinePossibleGroups()
{
    int alleleCounter = 0;
    for (Allele & allele : this->variant.getAlleles())
    {
        VariantRegions vRegions {
            createVariantRegions(allele)
        };
        determineVariantGroups(vRegions);
        if (allele.getName() != "REF")
        {
            this->variantAlleleNames[allele.getName()] = alleleCounter;
            ++alleleCounter;
        }
    } 
}

void VariantProfile::createVariantChromosomeStructures()
{
    this->chromosomeStructures.erase(this->chromosomeStructures.begin(), this->chromosomeStructures.end());
    for (Allele & allele : this->variant.getAlleles())
        createChromosomeStructures(allele);
}

void VariantProfile::initMasks()
{
    initReferenceMask();
    initVariantMask();
}

inline void VariantProfile::createIndexString(std::string & g, const std::unordered_set<int> & indexSet)
{
    std::vector<int> indices;
    for (const int & idx : indexSet)
        indices.push_back(idx);
    std::sort(indices.begin(), indices.end());
    g = "";
    for (const int & idx : indices)
        g += std::to_string(idx);
}

inline void VariantProfile::createIndexString(std::string & g, const std::unordered_set<std::string> & indexSet)
{
    std::vector<std::string> indices;
    for (const std::string & idx : indexSet)
        indices.push_back(idx);
    std::sort(indices.begin(), indices.end());
    g = "";
    for (const std::string & idx : indices)
        g += idx;
}

VariantRegions VariantProfile::createVariantRegions(Allele & allele)
{
    createChromosomeStructures(allele);
    std::unordered_map<std::string, JunctionRegion> chrStructures = this->chromosomeStructures[allele.getName()];

    VariantRegions vRegions;
    for (auto it : chrStructures)
    {
        const std::string & chr = it.first;
        JunctionRegion & jRegion = it.second;

        // group indices by distance
        std::vector<int32_t> indices;
        std::vector<int> localIndices;
        std::vector<char> type;

        for (uint32_t i = 0; i < jRegion.breakpointIndices.size(); ++i) 
        {
            indices.push_back(jRegion.breakpointIndices[i]);
            localIndices.push_back((int) i);
            type.push_back('B');
        }
        for (uint32_t i = 0; i < jRegion.junctionIndices.size(); ++i)
        {
            indices.push_back(jRegion.junctionIndices[i]);
            localIndices.push_back((int) i);
            type.push_back('J');
        }
        std::vector<int> helper(indices.size());
        std::iota(helper.begin(), helper.end(), 0);
        std::sort(helper.begin(), helper.end(), [&](const int & i, const int & j)
        {
            return (indices[i] < indices[j]);
        });

        std::sort(indices.begin(), indices.end());

        int regionIdx = 0;
        std::vector<std::vector<int>> grouping {std::vector<int> {0}};
        std::vector<int> regionDistances;
        for (uint32_t i = 1; i < indices.size(); ++i)
        {
            int32_t d = indices[i] - indices[i-1];
            if (type[i-1] == 'J')
                --d;

            if (d > (2*this->filterMargin))
            {
                ++regionIdx;
                grouping.push_back(std::vector<int>{(int) i});
                regionDistances.push_back(d - (2*this->filterMargin));
            } else {
                grouping[regionIdx].push_back(i);
            }
        }

        // go though all groups; define region, copy breakpoints, junctions 
        // and indices; store distances to other groups
        for (int g = 0; g <= regionIdx; ++g)
        {
            JunctionRegion newRegion;
            uint64_t minIdx {400000000}, maxIdx {0};
            char maxType{' '};
            for (auto & p : grouping[g])
            {
                uint64_t idx = indices[p];

                char & t = type[helper[p]];
                int & lI = localIndices[helper[p]];

                if (idx < minIdx)
                {
                    minIdx = idx;
                }
                if (idx > maxIdx)
                {
                    maxType = t;
                    maxIdx = idx;
                }

                if (t == 'B')
                {
                    newRegion.breakpoints.push_back(jRegion.breakpoints[lI]);
                    newRegion.breakpointIndices.push_back(jRegion.breakpointIndices[lI]);
                } 
                else if (t == 'J')
                {
                    newRegion.junctions.push_back(jRegion.junctions[lI]);
                    newRegion.junctionIndices.push_back(jRegion.junctionIndices[lI]);
                }
            }

            minIdx -= this->filterMargin;
            maxIdx += this->filterMargin;
            if (maxType == 'J') // need to make sure the correct side of junction is used
                ++maxIdx;
            newRegion.length = maxIdx - minIdx + 1;
            
            // store the distance to the other groups
            int32_t nextDistance{-1}, lastDistance{-1};
            if (g > 0)
                lastDistance = regionDistances[g-1];
            if (g < regionIdx)
                nextDistance = regionDistances[g];
            vRegions.distanceToNext.push_back(nextDistance);
            vRegions.distanceFromLast.push_back(lastDistance);
            newRegion.chromosome = chr;

            // define the region(s)
            newRegion.regions = getRegionsFromIndices(jRegion, minIdx, maxIdx);

            // adjust the breakpoint and junction indices w.r.t. region start
            for (auto & j : newRegion.junctionIndices)
                j -= minIdx;
            for (auto & b: newRegion.breakpointIndices)
                b -= minIdx;
            vRegions.regions.push_back(newRegion);
        }
    }
    
    return vRegions;
}

void VariantProfile::createChromosomeStructures(Allele & allele)
{
    // find unique chromosomes
    std::vector<Breakpoint> & breakpoints = this->variant.getAllBreakpoints();
    std::vector<Junction> & junctions = allele.getNovelJunctions();

    std::unordered_set<std::string> chromosomes;
    for (auto & bp : breakpoints)
        chromosomes.insert(bp.getReferenceName());

    // go over all involved chromosomes
    for (auto & chr : chromosomes)
    {
        JunctionRegion chrStructure;
        // determine breakpoints and junctions on the chromosome
        std::vector<int> bpIndices;
        std::vector<int> njIndices;
        int32_t regionLength = 0;
        for (uint32_t i = 0; i < breakpoints.size(); ++i)
            if (breakpoints[i].getReferenceName() == chr)
                bpIndices.push_back((int) i);
        for (uint32_t i = 0; i < junctions.size(); ++i)
            if (junctions[i].getVariantRefName() == chr)
                njIndices.push_back((int) i);
        std::sort(bpIndices.begin(), bpIndices.end(), [&](const uint32_t & i, const uint32_t & j)
        {
            return breakpoints[i].getPosition() < breakpoints[j].getPosition();
        });

        if (njIndices.size() > 0)
        {
            // create a set of regions describing the chromosome structure
            // left
            bool isReverse {junctions[njIndices[0]].getDirectionLeft() > 0};
            chrStructure.regions.push_back(
                GenomicRegion {
                    junctions[njIndices[0]].getRefNameLeft(),
                    (isReverse ? junctions[njIndices[0]].getPositionLeft() : 0),
                    (isReverse ? this->cInfo.getContigLengths()[junctions[njIndices[0]].getRefNameLeft()] : junctions[njIndices[0]].getPositionLeft()),
                    isReverse
                }
            );
            regionLength = chrStructure.regions[0].getRegionEnd() - chrStructure.regions[0].getRegionStart();
            chrStructure.junctionIndices.push_back(regionLength);
            chrStructure.length = regionLength + 1;
            chrStructure.junctions.push_back(junctions[njIndices[0]]);
            // middle
            for (uint32_t j = 0; j + 1 < njIndices.size(); ++j)
            {
                if (junctions[njIndices[j]].getRefNameRight() != junctions[njIndices[j + 1]].getRefNameLeft())
                    throw std::runtime_error("ERROR: RefNames of adjacent junctions do not match!");
                if (junctions[njIndices[j]].getDirectionRight() == junctions[njIndices[j + 1]].getDirectionLeft())
                    throw std::runtime_error("ERROR: Directions of adjacent junctions are incompatible!");
                isReverse = (junctions[njIndices[j]].getDirectionRight() < 0);
                chrStructure.regions.push_back(
                    GenomicRegion {
                        junctions[njIndices[j]].getRefNameRight(),
                        (isReverse ? junctions[njIndices[j+1]].getPositionLeft() : junctions[njIndices[j]].getPositionRight()),
                        (isReverse ? junctions[njIndices[j]].getPositionRight() : junctions[njIndices[j+1]].getPositionLeft()),
                        isReverse
                    }
                );
                regionLength = chrStructure.regions[j+1].getRegionEnd() - chrStructure.regions[j+1].getRegionStart() + 1;
                chrStructure.junctionIndices.push_back(
                    regionLength + chrStructure.junctionIndices[j]
                );
                chrStructure.length += regionLength;
                chrStructure.junctions.push_back(junctions[njIndices[j+1]]);
            }
            // right
            int tempIdx {njIndices[njIndices.size() - 1]};
            isReverse = junctions[tempIdx].getDirectionRight() < 0;

            int32_t cLen = 0;
            if (this->cInfo.getContigLengths().find(junctions[tempIdx].getRefNameRight()) == this->cInfo.getContigLengths().end()) {
                cLen = 300000000;
                std::cerr << "Error: " << junctions[tempIdx].getRefNameRight() << " not found in contigs" << std::endl;
            } else {
                cLen = this->cInfo.getContigLengths()[junctions[tempIdx].getRefNameRight()];
            }
            chrStructure.regions.push_back(
                GenomicRegion {
                    junctions[tempIdx].getRefNameRight(),
                    (isReverse ? 0 : junctions[tempIdx].getPositionRight()),
                    (isReverse ? junctions[tempIdx].getPositionRight() : cLen),
                    isReverse
                }
            );
            regionLength = chrStructure.regions[njIndices.size()].getRegionEnd() - chrStructure.regions[njIndices.size()].getRegionStart() + 1;
            chrStructure.length += regionLength;
        } 
        else 
        {
            chrStructure = JunctionRegion {
                std::vector<GenomicRegion> {
                    GenomicRegion{
                        chr,
                        0,
                        cInfo.getContigLengths()[chr],
                        false
                    }
                },
                std::vector<int>(),
                std::vector<int>(),
                std::vector<Junction>(),
                std::vector<Breakpoint>(),
                cInfo.getContigLengths()[chr]
            };
        }

        // Place all breakpoints and assign corresponding indices
        for (auto & b : bpIndices)
            insertBreakpoint(breakpoints[b], chrStructure);
        
        this->chromosomeStructures[allele.getName()][chr] = chrStructure;
    }
}

void VariantProfile::insertBreakpoint(Breakpoint & bp, JunctionRegion & jRegion)
{
    // find fitting regions
    int currentIdx = 0;
    for (uint32_t j = 0; j < jRegion.regions.size(); ++j)
    {
        // adjust region to not include breakpoint locations
        GenomicRegion tempRegion {
            jRegion.regions[j].getReferenceName(),
            jRegion.regions[j].getRegionStart() + 1,
            jRegion.regions[j].getRegionEnd() - 1,
            jRegion.regions[j].isReverse()
        };
        if (tempRegion.overlaps(GenomicRegion {
            bp.getReferenceName(),
            bp.getPosition(),
            bp.getPosition()
            })
        )
        {
            // determine exact index and add breakpoint and index to struct
            int d;
            if (jRegion.regions[j].isReverse())
                d = jRegion.regions[j].getRegionEnd() - bp.getPosition();
            else
                d = bp.getPosition() - jRegion.regions[j].getRegionStart();
            int index = currentIdx + d;
            jRegion.breakpoints.push_back(bp);
            jRegion.breakpointIndices.push_back(index);
        }
        currentIdx += jRegion.regions[j].getRegionEnd() - jRegion.regions[j].getRegionStart() + 1;
    }
}

std::vector<GenomicRegion> VariantProfile::getRegionsFromIndices(JunctionRegion & jRegion, int minIdx, int maxIdx)
{
    std::vector<GenomicRegion> regions;
    int rStartIdx {-1}, rEndIdx{-1};
    for (auto & r : jRegion.regions)
    {
        rStartIdx = rEndIdx + 1;
        rEndIdx += r.getRegionEnd() - r.getRegionStart() + 1;
        
        if (rEndIdx > minIdx && rStartIdx < maxIdx)
        {
            // first partial region
            if (minIdx >= rStartIdx && minIdx <= rEndIdx && maxIdx > rEndIdx)
            {
                GenomicRegion temp = r;
                int d = minIdx - rStartIdx;
                if (r.isReverse())
                    temp.setRegionEnd(r.getRegionEnd() - d);
                else 
                    temp.setRegionStart(r.getRegionStart() + d);
                regions.push_back(temp);
            }

            // required complete regions
            if (rEndIdx < maxIdx && rStartIdx > minIdx)
                regions.push_back(r);

            // middle partial region
            if (minIdx > rStartIdx && maxIdx < rEndIdx)
            {
                GenomicRegion temp = r;
                int dMin = minIdx - rStartIdx;
                int dMax = rEndIdx - maxIdx;
                if (r.isReverse())
                {
                    temp.setRegionStart(r.getRegionStart() + dMax);
                    temp.setRegionEnd(r.getRegionEnd() - dMin);
                } 
                else 
                {
                    temp.setRegionStart(r.getRegionStart() + dMin);
                    temp.setRegionEnd(r.getRegionEnd() - dMax);
                }
                regions.push_back(temp);
            }

            // last partial region
            if (maxIdx <= rEndIdx && maxIdx >= rStartIdx && minIdx < rStartIdx)
            {
                GenomicRegion temp = r;
                int d = maxIdx - rStartIdx;
                if (r.isReverse())
                    temp.setRegionStart(r.getRegionEnd() - d);
                else
                    temp.setRegionEnd(r.getRegionStart() + d);
                regions.push_back(temp);
            }
        }
    }
    return regions;
}

void VariantProfile::determineVariantGroups(VariantRegions & vRegions)
{
    std::unordered_set<std::string> groups;
    findPairAttributes(groups, vRegions);
    determineSpanningGroups(groups, vRegions);
    determineSplitGroups(groups, vRegions);

    int counter = this->variantGroups.size();
    for (auto & g : groups)
    {
        if (this->variantGroups.find(g) == this->variantGroups.end())
        {
            this->variantGroups[g] = counter;
            this->groupOccurs[g] = false;
            ++counter;
        }
    }
}

void VariantProfile::findPairAttributes(std::unordered_set<std::string> & groups, VariantRegions & vRegions)
{
    groups.insert("FR");
    groups.insert("RR");
    groups.insert("FF");
    
    // go over all combinations of regions
    for (int r1Idx = 0; r1Idx < (int) vRegions.regions.size(); ++r1Idx)
    {
        auto & jRegion = vRegions.regions[r1Idx];
        int32_t d = 0;

        for (int r2Idx = r1Idx; r2Idx >= 0; --r2Idx)
        {
            // only consider combination if distance between regions is not larger than max insert size
            if (r2Idx < r1Idx)
                d += vRegions.distanceToNext[r2Idx];
            if (d > this->sMax)
                break;

            auto & kRegion = vRegions.regions[r2Idx];
            // all combinations of junctions in the regions
            for (int j = 0; j < (int) jRegion.junctions.size(); ++j)
            {
                for (int k = (int) kRegion.junctions.size() - 1; k >= 0; --k)
                {
                    if (r1Idx == r2Idx && k > j)
                        continue;
                    // check distance of junctions
                    if (r1Idx != r2Idx && (kRegion.length - 1 - kRegion.junctionIndices[k]) + d + jRegion.junctionIndices[j] > this->sMax)
                        break;
                    if (r1Idx == r2Idx && (jRegion.junctionIndices[j] - kRegion.junctionIndices[k]) > this->sMax)
                        break;

                    std::vector<Breakpoint> kBps = {kRegion.junctions[k].leftSideToBreakpoint(0), kRegion.junctions[k].rightSideToBreakpoint(1)};
                    std::vector<Breakpoint> jBps = {jRegion.junctions[j].leftSideToBreakpoint(0), jRegion.junctions[j].rightSideToBreakpoint(1)};
                    for (auto & kB : kBps)
                    {
                        for (auto & jB : jBps)
                        {
                            if (cInfo.globalPositions.find(jB.getReferenceName()) == cInfo.globalPositions.end())
                                std::cerr << "Could not find contig " << jB.getReferenceName() << " in contigInfo." << std::endl;
                            if (cInfo.globalPositions.find(kB.getReferenceName()) == cInfo.globalPositions.end())
                                std::cerr << "Could not find contig " << kB.getReferenceName() << " in contigInfo." << std::endl;
                            
                            int64_t s = (int64_t) this->cInfo.globalPositions[jB.getReferenceName()] + jB.getPosition() -
                                ((int64_t) this->cInfo.globalPositions[kB.getReferenceName()] + kB.getPosition());

                            this->sMinMapped = std::min(this->sMinMapped, s);
                            this->sMaxMapped = std::max(this->sMaxMapped, std::abs(s));
                        }
                    }
                }
            }
        }
    }
    this->sMinMapped -= this->sMax;
    this->sMaxMapped += this->sMax;
}


bool VariantProfile::isPossible(std::vector<int> positions)
{
    if (positions.size() == 0)
        return false;

    std::sort(positions.begin(), positions.end());
    if (positions[positions.size() - 1] > this->sMax)
        return false;

    uint32_t j = 0;
    for (; j < positions.size(); ++j)
        if (positions[j] > this->readLength)
            break;
    int k = (int) positions.size() - 1;
    for (; k >= 0; --k)
        if (positions[positions.size() - 1] - positions[k] > this->readLength || k <= (int) j + 1)
            break;
    
    return (k <= (int) j + 1);
}

std::vector<std::vector<int>> VariantProfile::getSubsets(std::vector<int> positions, std::vector<int> ids, uint32_t idx)
{
    std::vector<std::vector<int>> subsets;
    if (positions.size() == 0)
	    return subsets;

    for (uint32_t i = idx; i < positions.size(); ++i)
    {
        std::vector<int> tempPositions = positions;
        std::vector<int> tempIDs = ids;
        tempPositions.erase(tempPositions.begin() + i);
        tempIDs.erase(tempIDs.begin() + i);
        std::vector<std::vector<int>> tempSets = getSubsets(tempPositions, tempIDs, i);
        for (auto & s : tempSets)
            subsets.push_back(s);
    }

    if (isPossible(positions))
        subsets.push_back(ids);
    return subsets;
}

void VariantProfile::determineSplitGroups(std::unordered_set<std::string> & groups, VariantRegions & vRegions)
{
    std::vector<std::unordered_set<int>> indexGroups;

    // consider all regions
    for (uint32_t v = 0; v < vRegions.regions.size(); ++v)
    {
        JunctionRegion & jRegion = vRegions.regions[v];

        // consider all junctions as possible starting positions for read pairs
        for (uint32_t i = 0; i < jRegion.junctions.size(); ++i)
        {
            std::vector<int> ids;
            std::vector<int32_t> positions;
            bool extend = true;
            ids.push_back(jRegion.junctions[i].getID());
            positions.push_back(0);
            uint32_t j = i + 1;
            while (j < jRegion.junctions.size())
            {
                if (jRegion.junctionIndices[j] - jRegion.junctionIndices[i] > this->sMax)
                {
                    extend = false;
                    break;
                }
                ids.push_back(jRegion.junctions[j].getID());
                positions.push_back(jRegion.junctionIndices[j] - jRegion.junctionIndices[i]);
		        ++j;
            }
            int32_t d = jRegion.length - 1 - jRegion.junctionIndices[i];
            uint32_t w = v + 1;
            while (extend && vRegions.distanceToNext[w - 1] >= 0)
            {
                d += vRegions.distanceToNext[w - 1];
                for (uint32_t k = 0; k < vRegions.regions[w].junctionIndices.size(); ++k)
                {
                    if (d + vRegions.regions[w].junctionIndices[k] > this->sMax)
                    {
                        extend = false;
                        break;
                    }
                    ids.push_back(vRegions.regions[w].junctions[k].getID());
                    positions.push_back(d + vRegions.regions[w].junctionIndices[k]);
                }
                d += vRegions.regions[w].length - 1;
                ++w;
            }
            // collect all possible split read groups
            std::vector<std::vector<int>> indexSubsets = getSubsets(positions, ids, 0);
            for (auto & s : indexSubsets)
            {
                std::unordered_set<int> tempSet;
                for (auto & idx : s)
                    tempSet.insert(idx);
                indexGroups.push_back(tempSet);
            }
        }
    }

    // create final unique group strings from junction indices
    std::string splitString;
    for (uint32_t i = 1; i <= indexGroups.size(); ++i)
    {
        createIndexString(splitString, indexGroups[i-1]);
        splitString = "split_" + splitString;
        groups.insert(splitString);
    }
    groups.insert("split_ambiguous");
}


void VariantProfile::determineSpanningGroups(std::unordered_set<std::string> & groups, VariantRegions & vRegions)
{
    std::vector<std::unordered_set<int>> indexGroups;

    // consider all regions
    for (uint32_t v = 0; v < vRegions.regions.size(); ++v)
    {
        JunctionRegion & bRegion = vRegions.regions[v];
        // consider all junctions as possible starting positions for read pairs
        for (uint32_t i = 0; i < bRegion.breakpoints.size(); ++i)
        {
            std::vector<int> ids;
            std::vector<int32_t> positions;
            bool extend = true;
            ids.push_back(bRegion.breakpoints[i].getID());
            positions.push_back(0);
            uint32_t j = i + 1;
            while (j < bRegion.breakpoints.size())
            {
                if (bRegion.breakpointIndices[j] - bRegion.breakpointIndices[i] > this->sMax)
                {
                    extend = false;
                    break;
                }
                ids.push_back(bRegion.breakpoints[j].getID());
                positions.push_back(bRegion.breakpointIndices[j] - bRegion.breakpointIndices[i]);
		        ++j;
            }
            int32_t d = bRegion.length - 1 - bRegion.breakpointIndices[i];

            int w = v + 1;
            while (extend && vRegions.distanceToNext[w - 1] >= 0)
            {
                d += vRegions.distanceToNext[w - 1];
                for (uint32_t k = 0; k < vRegions.regions[w].breakpointIndices.size(); ++k)
                {
                    if (d + vRegions.regions[w].breakpointIndices[k] > this->sMax)
                    {
                        extend = false;
                        break;
                    }
                    ids.push_back(vRegions.regions[w].breakpoints[k].getID());
                    positions.push_back(d + vRegions.regions[w].breakpointIndices[k]);
                }
                d += vRegions.regions[w].length - 1;
                w++;
            }
            
            // collect all possible spanning read groups
            std::vector<std::vector<int>> indexSubsets = getSubsets(positions, ids, 0);
            for (auto & s : indexSubsets)
            {
                std::unordered_set<int> tempSet;
                for (auto & idx : s)
                    tempSet.insert(idx);
                indexGroups.push_back(tempSet);
            }
        }
    }

    // create final unique group strings from junction indices
    std::string spanningString;
    std::string bridgeString;
    for (uint32_t i = 1; i <= indexGroups.size(); ++i)
    {
        createIndexString(spanningString, indexGroups[i-1]);
        bridgeString  = "bridging_" + spanningString;
        spanningString = "spanning_" + spanningString;
        groups.insert(spanningString);
        groups.insert(bridgeString);
    }
}


void VariantProfile::initReferenceMask()
{ 
    this->referenceMask = Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t>(this->variantGroups.size(), this->sMaxMapped - this->sMinMapped + 1);
    this->referenceMask.reserve(Eigen::VectorXi::Constant(this->variantGroups.size(), this->sMax - this->sMin + 1)); 
}

void VariantProfile::initVariantMask()
{ 
    for (uint32_t i = 0; i + 1 < this->variant.getAlleles().size(); ++i)
    {
        this->variantMask.push_back(
            std::vector<Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t>>(
                this->sMax - this->sMin + 1, 
                Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t>(this->variantGroups.size(), this->sMaxMapped - this->sMinMapped + 1)
            )
        );
        for (int j = 0; j < this->sMax - this->sMin + 1; ++j)
            this->variantMask[i][j].reserve(Eigen::VectorXi::Constant(this->variantGroups.size(), (this->variant.getAllJunctions().size() + 1) *(this->sMax - this->sMin + 1)));
    } 
}


void VariantProfile::calculateAlleleMasks()
{ 
    int status = -1;
    
    for (Allele & allele : this->variant.getAlleles())
    {
        for (std::string cName : allele.getChromosomeNames())
        {
            VariantMapManager & chrMaps = allele.getChromosomeMap(cName);

            std::vector<VariantMap> maps = chrMaps.getMaps();

            for (uint32_t j = 0; j < maps.size(); ++j)
            {
                VariantMap & map = maps[j];
                
                int mapLength = map.totalLength;
                for (int i = 0; i < mapLength; ++i)
                {
                    for (int64_t s = this->sMax; s >= this->sMin; --s)
                    {
                        if (i+s >= mapLength)
                            continue;
                        addSimulatedTemplateToMask(
                            status, map, i, s, allele
                        );
                        if (status < 0)
                            break;
                    }
                }
            }
        }
    }

    this->referenceMask.makeCompressed();
    for (uint32_t i = 0; i + 1 < this->variant.getAlleles().size(); ++i)
        for (int j = 0; j < this->sMax - this->sMin + 1; ++j)
            this->variantMask[i][j].makeCompressed(); 
}


inline void VariantProfile::addSimulatedTemplateToMask(int & status, VariantMap & map, int i, int64_t s, Allele & allele)
{
    ReadTemplate sT = map.simulateTemplate(i, s, this->readLength, this->cInfo);
    if (!sT.isInterestingReadPair(this->filter))
    {
        status = -1;
        return;
    }

    sT.findSpanningReads(this->variant.getAllBreakpoints());
    sT.findBridgedBreakpoints(this->variant.getAllBreakpoints());
    if (sT.containsSuspectedSplit())
        sT.findSplitReads(allele.getNovelJunctions(), this->chromosomeStructures, this->sMax);
    sT.determineLocationStrings();

    std::string orientation = sT.getOrientation();
    int64_t insertSize = sT.getInsertSize();
    std::string junctionString = sT.getJunctionString();
    std::string breakpointString = sT.getBpSpanString();
    std::string bridgeString = sT.getBpBridgeString();

    status = 1;
    addValueToMask(
        allele, 
        s, 
        insertSize, 
        orientation, 
        junctionString, 
        breakpointString, 
        bridgeString
        );

    return;
}

inline void VariantProfile::addValueToMask(Allele & allele, int64_t sOld, int64_t sNew, std::string & orientation, std::string & jString, std::string & bpString, std::string & bridgeString)
{
    std::string group = GenotypeDistribution::determineGroup(orientation, jString, bpString, bridgeString);

    if (this->variantGroups.find(group) == this->variantGroups.end()) {
        std::cout << "Warning: Could not assign read pair to any group. Adding group." << std::endl;
        std::cout << "Variant: " << this->variant.getName() << std::endl; 
        std::cout << "Group: " << group << std::endl;

        // 
        // in case of group error: conservative matrix resize; try to fix the underlying issue (determineSplitGroups) and 
        // allocate memory a little bit more generously
        //
        addGroupToMasks(group);
    }

    if (sNew < this->sMinMapped || sNew >= this->sMaxMapped)
	    std::cerr << "Insert size error. Profile boundaries: " << this->sMinMapped << "\t" << this->sMaxMapped << "\tMapped insert size: " << sNew << std::endl;

    int gIdx = this->variantGroups[group];
    this->groupOccurs[group] = true;

    if (allele.getName() == "REF") 
        this->referenceMask.coeffRef(gIdx, sOld-this->sMinMapped) += 1;
    else 
        this->variantMask[this->variantAlleleNames[allele.getName()]][sOld-this->sMin].coeffRef(gIdx, sNew - this->sMinMapped) += 1;
 
    return;
}

void VariantProfile::addGroupToMasks(std::string group)
{
    if (this->variantGroups.find(group) != this->variantGroups.end())
        return;
    int nG = this->variantGroups.size();
    this->variantGroups[group] = nG;
    this->groupOccurs[group] = true;

    this->referenceMask.conservativeResize(this->variantGroups.size(), this->referenceMask.cols());
    this->referenceMask.reserve(Eigen::VectorXi::Constant(this->variantGroups.size(), this->sMax - this->sMin + 1));
    
    for (uint32_t i = 0; i < this->variantMask.size(); ++i)
    {
        for (uint32_t j = 0; j < this->variantMask[i].size(); ++j)
        {
            this->variantMask[i][j].conservativeResize(this->variantGroups.size(), this->variantMask[i][j].cols());
            this->variantMask[i][j].reserve(Eigen::VectorXi::Constant(this->variantGroups.size(), (this->variant.getAllJunctions().size() + 1) *(this->sMax - this->sMin + 1)));
        }
    }
}

const Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t> & VariantProfile::getVariantMask(int s)
{
    return this->variantMask[0][s];
}

const Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t> & VariantProfile::getVariantMask(const std::string & alleleName, int s)
{
    return this->variantMask[this->variantAlleleNames[alleleName]][s];
}

const Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t> & VariantProfile::getReferenceMask()
{
    return this->referenceMask;
}

void VariantProfile::calculateGenotypeDistributions(std::unordered_map<std::string, GenotypeDistribution> & distributions, LibraryDistribution & libraryDistribution, float eps)
{   
    std::vector<std::string> alleleNames(this->variantAlleleNames.size() + 1, "");
    alleleNames[0] = "REF";
    for (auto & vN : this->variantAlleleNames)
        alleleNames[vN.second + 1] = vN.first;

    // create allele distributions
    std::vector<Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t>> variantAlleleDists(
        alleleNames.size(), 
        Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t>(this->variantGroups.size(), this->sMaxMapped - this->sMinMapped + 1)
    );

    // reserve space for allele matrices
    for (uint32_t i = 0; i < alleleNames.size(); ++i)
        variantAlleleDists[i].reserve(Eigen::VectorXi::Constant(this->variantGroups.size(), (this->variant.getAllJunctions().size() + 1) * (this->sMax - this->sMin + 1)));

    // fill matrices
    for (int64_t s = this->sMin; s <= this->sMax; ++s)
    {
        float p = libraryDistribution.getProbability((int) s);
        for (uint32_t v = 0; v < alleleNames.size(); ++v)
        {	
	        if (alleleNames[v] == "REF")
                for (int i = 0; i < variantAlleleDists[0].rows(); ++i)
                    variantAlleleDists[0].insert(i, s - this->sMinMapped) = this->referenceMask.coeff(i, s-this->sMinMapped) * p;
	        else
                variantAlleleDists[v] += this->variantMask[this->variantAlleleNames[alleleNames[v]]][s - this->sMin] * p; 
	    }
    }

    // mix distributions
    float majorFactor = 1. - eps;
    float minorFactor = eps / variantAlleleDists.size();

    std::vector<Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t>> variantGtDists;
    std::vector<std::string> gtNames;

    Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t> temp(variantAlleleDists[0].rows(), variantAlleleDists[0].cols());
    for (uint32_t i = 0; i < variantAlleleDists.size(); ++i)
    {
        temp = (2 * variantAlleleDists[i] * majorFactor).pruned();
        for (uint32_t j = 0; j < variantAlleleDists.size(); ++j)
        {
            if (j == i)
                continue;
            temp += (minorFactor * variantAlleleDists[j]).pruned();
        }
        
        variantGtDists.push_back(temp.pruned());
        gtNames.push_back(alleleNames[i] + "/" + alleleNames[i]);

        for (uint32_t j = i + 1; j < variantAlleleDists.size(); ++j)
        {
            temp  = ((variantAlleleDists[i] + variantAlleleDists[j]) * majorFactor).pruned();
            for (uint32_t k = 0; k < variantAlleleDists.size(); ++k)
            {
                if (k == i || k == j)
                    continue;
                temp += (minorFactor * variantAlleleDists[k]).pruned();
            }
            
            variantGtDists.push_back(temp.pruned());
            gtNames.push_back(alleleNames[i] + "/" + alleleNames[j]);
        }
    }

    // scale distributions
    for (auto & dist : variantGtDists)
        dist /= dist.sum();

    // create GenotypeDistributions and map
    float minProb = 1;
    for (uint32_t i = 0; i < variantGtDists.size(); ++i) 
    {
        distributions[gtNames[i]] = GenotypeDistribution(variantGtDists[i], this->variantGroups, this->sMinMapped, this->sMaxMapped);
        if (distributions[gtNames[i]].getMinProbability() > 0)
            minProb = std::min(minProb, distributions[gtNames[i]].getMinProbability());
    }
    for (uint32_t i = 0; i < variantGtDists.size(); ++i)
        distributions[gtNames[i]].setMinProbability(minProb);
    
    return;
}

ReadPairFilter & VariantProfile::getFilter()
{
    return this->filter;
}

complexVariant & VariantProfile::getVariant()
{
    return this->variant;
}

const ContigInfo & VariantProfile::getContigInfo()
{
    return this->cInfo;
}

std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> & VariantProfile::getChromosomeStructures()
{
    return this->chromosomeStructures;
}

void VariantProfile::writeProfile(std::string filename)
{
    std::ofstream stream(filename, std::ios_base::binary | std::ios_base::out);
    if (!stream.is_open())
    {
        std::string msg = "Could not open file " + filename + " to write Profile";
        throw std::runtime_error(msg.c_str());
    }

    // write magic string
    stream.write("GENOTYPER\1", 10);
    
    // write variant name
    int vnSize = this->variant.getName().size() + 1;
    stream.write(reinterpret_cast<const char *>(&vnSize), sizeof(int));
    stream.write(this->variant.getName().c_str(), vnSize - 1);
    stream.write("\0", 1);

    // write variant file name
    std::string vName = this->variant.getVariantFileName();
    int vfSize = vName.size() + 1;
    stream.write(reinterpret_cast<const char *>(&vfSize), sizeof(int));
    stream.write(vName.c_str(), vfSize - 1);
    stream.write("\0", 1);
    
    // write read length, sMin and sMax assumed for library distibution
    stream.write(reinterpret_cast<const char *>(&this->readLength), sizeof(int));
    stream.write(reinterpret_cast<const char *>(&this->sMin), sizeof(int));
    stream.write(reinterpret_cast<const char *>(&this->sMax), sizeof(int));

    // write min and max insert size in profiles
    stream.write(reinterpret_cast<const char *>(&this->sMinMapped), sizeof(int64_t));
    stream.write(reinterpret_cast<const char *>(&this->sMaxMapped), sizeof(int64_t));

    // write filter margin
    stream.write(reinterpret_cast<const char *>(&this->filterMargin), sizeof(int));

    // write chromosome names and lengths
    // write number of chromosomes
    uint32_t nChrom = this->cInfo.cNames.size();
    stream.write(reinterpret_cast<const char *>(&nChrom), sizeof nChrom);
    for (uint32_t i = 0; i < this->cInfo.cNames.size(); ++i)
    {
        uint32_t cNameLen = this->cInfo.cNames[i].size() + 1;
        stream.write(reinterpret_cast<const char *>(&cNameLen), sizeof cNameLen);
        stream.write(this->cInfo.cNames[i].c_str(), cNameLen - 1);
        stream.write("\0", 1);
        stream.write(reinterpret_cast<const char *>(&this->cInfo.cLengths[i]), sizeof(int32_t));
    }

    // write number of read pair groups and their names
    int nG = this->variantGroups.size();
    stream.write(reinterpret_cast<const char *>(&nG), sizeof(int));
    
    std::vector<std::string> colnames(nG, "");
    for (int i = 0; i < nG; )
    {
        for (auto & g : this->variantGroups)
        {
            if (g.second == i)
            {
                colnames[i] = g.first;
                ++i;
            }
        }
    }
    for (std::string & g : colnames)
    {
        int size = g.size() + 1;
        stream.write(reinterpret_cast<const char *>(&size), sizeof(int));
        stream.write(g.c_str(), size - 1);
        stream.write("\0", 1);
    }

    // write reference allele name
    stream.write("REF\0", 4);

    // write reference matrix
    // this->referenceMask.makeCompressed()  ;
    int64_t m      = referenceMask.rows()     ;
    int64_t n      = referenceMask.cols()     ;
    int64_t nnzS   = referenceMask.nonZeros() ;
    int64_t outerS = referenceMask.outerSize();
    int64_t innerS = referenceMask.innerSize();
    stream.write(reinterpret_cast<const char *>(&m), sizeof(int64_t));
    stream.write(reinterpret_cast<const char *>(&n), sizeof(int64_t));
    stream.write(reinterpret_cast<const char *>(&nnzS), sizeof(int64_t));
    stream.write(reinterpret_cast<const char *>(&outerS), sizeof(int64_t));
    stream.write(reinterpret_cast<const char *>(&innerS), sizeof(int64_t));

    stream.write(reinterpret_cast<const char *>(referenceMask.valuePtr()), sizeof(float)*nnzS);
    stream.write(reinterpret_cast<const char *>(referenceMask.outerIndexPtr()), sizeof(int64_t)*outerS);
    stream.write(reinterpret_cast<const char *>(referenceMask.innerIndexPtr()), sizeof(int64_t)*nnzS);

    // write number of variant alleles and their names
    int nA = this->variantAlleleNames.size();
    stream.write(reinterpret_cast<const char *>(&nA), sizeof(int));
    std::vector<std::string> alleleNames;
    for (int i = 0; i < nA; )
    {
        for (auto & a : this->variantAlleleNames)
        {
            if (a.second == i)
            {
                alleleNames.push_back(a.first);
                ++i;
            }
        }
    }
    for (std::string & a : alleleNames)
    {
        int size = a.size() + 1;
        stream.write(reinterpret_cast<const char *>(&size), sizeof(int));
        stream.write(a.c_str(), size - 1);
        stream.write("\0", 1);
    }

    // write variant matrices
    for (uint32_t i = 0; i < this->variantMask.size(); ++i)
    {
        for (uint32_t s = 0; s < this->variantMask[i].size(); ++s)
        {
            // this->variantMask[i][s].makeCompressed();

            int64_t m      = this->variantMask[i][s].rows()     ;
            int64_t n      = this->variantMask[i][s].cols()     ;
            int64_t nnzS   = this->variantMask[i][s].nonZeros() ;
            int64_t outerS = this->variantMask[i][s].outerSize();
            int64_t innerS = this->variantMask[i][s].innerSize();

            stream.write(reinterpret_cast<const char *>(&m), sizeof(int64_t))     ;
            stream.write(reinterpret_cast<const char *>(&n), sizeof(int64_t))     ;
            stream.write(reinterpret_cast<const char *>(&nnzS), sizeof(int64_t))  ;
            stream.write(reinterpret_cast<const char *>(&outerS), sizeof(int64_t));
            stream.write(reinterpret_cast<const char *>(&innerS), sizeof(int64_t));

            stream.write(reinterpret_cast<const char *>(this->variantMask[i][s].valuePtr()), sizeof(float)*nnzS)     ;
            stream.write(reinterpret_cast<const char *>(this->variantMask[i][s].outerIndexPtr()), sizeof(int64_t)*outerS);
            stream.write(reinterpret_cast<const char *>(this->variantMask[i][s].innerIndexPtr()), sizeof(int64_t)*nnzS)  ;
        }
    }
    
    stream.close();
}

void VariantProfile::readProfile(std::string filename)
{
    // open file
    std::ifstream stream(filename, std::ios_base::binary | std::ios_base::in);
    if (!stream.is_open())
    {
        std::string msg = "Could not open profile " + filename + " for reading";
        throw std::runtime_error(msg.c_str());
    }

    // compare magic string
    char * tempString = new char[10];
    stream.read(tempString, 10);
    if (tempString[9] != '\1')
        throw std::runtime_error("Magic string does not match. Wrong / corrupted profile?");
    tempString[9] = '\0';
    if (std::strcmp(tempString, "GENOTYPER") != 0)
        throw std::runtime_error("Magic string does not match. Wrong / corrupted profile?");
    delete[] tempString;

    // get variant name
    int l;
    stream.read(reinterpret_cast<char *>(&l), sizeof(int));
    tempString = new char[l];
    stream.read(tempString, l);
    if (tempString[l-1] != '\0')
        throw std::runtime_error("");
    this->name = std::string(tempString);
    delete[] tempString;
    
    // get variant path
    stream.read(reinterpret_cast<char *>(&l), sizeof(int));
    tempString = new char[l];
    stream.read(tempString, l);
    if (tempString[l-1] != '\0')
        throw std::runtime_error("");
    std::string variantPath(tempString);
    delete[] tempString;    

    // get read length, insert minima and maxima and filter margin
    stream.read(reinterpret_cast<char *>(&this->readLength), sizeof(int));
    stream.read(reinterpret_cast<char *>(&this->sMin), sizeof(int));
    stream.read(reinterpret_cast<char *>(&this->sMax), sizeof(int));
    stream.read(reinterpret_cast<char *>(&this->sMinMapped), sizeof(int64_t));
    stream.read(reinterpret_cast<char *>(&this->sMaxMapped), sizeof(int64_t));
    stream.read(reinterpret_cast<char *>(&this->filterMargin), sizeof(int));

    // chromosome names and lengths
    this->cInfo = ContigInfo();
    uint32_t nChrom;
    stream.read(reinterpret_cast<char * >(&nChrom), sizeof(uint32_t));
    for (uint32_t i = 0; i < nChrom; ++i)
    {
        uint32_t l;
        int cLength;

        stream.read(reinterpret_cast<char *>(&l), sizeof(uint32_t));

        char * tempString = new char[l];
        stream.read(reinterpret_cast<char *>(tempString), l);
        std::string cName = std::string(tempString);
        delete[] tempString;
        
        stream.read(reinterpret_cast<char *>(&cLength), sizeof(int32_t));
        this->cInfo.cNames.push_back(cName);
        this->cInfo.cLengths.push_back(cLength);
    }
    this->cInfo.calculateGlobalContigPositions();

    // try load the actual variant
    this->variantPresent = loadVariantStructure(variantPath, this->name);
    //
    // Maybe I should get rid of this dependency and write the filter to the profile?
    //
    if (variantPresent)
    {
        this->variant.setFilterMargin(this->filterMargin);
        this->filter = ReadPairFilter(this->variant.getAllBreakpoints(), this->filterMargin, 0);
    } 
    else 
    {
        std::string error = "Variant file not found in location given in the profile. Abort.\n";
        throw std::runtime_error(error.c_str());
    }

    // get read pair groups
    int nG;
    stream.read(reinterpret_cast<char *>(&nG), sizeof(int));
    std::unordered_map<std::string, int>().swap(this->variantGroups);
    for (int i = 0; i < nG; ++i)
    {
        stream.read(reinterpret_cast<char *>(&l), sizeof(int));
        tempString = new char[l];
        stream.read(tempString, l);
        if (tempString[l-1] != '\0')
            throw std::runtime_error("");
        std::string g(tempString);
        this->variantGroups[g] = i;
        delete[] tempString;
    }

    // get reference allele name
    tempString = new char[4];
    stream.read(tempString, 4);
    if (tempString[3] != '\0')
        throw std::runtime_error("");
    std::string refName(tempString);
    delete [] tempString;
    
    // get reference mask
    int64_t m, n, nnzS, innerS, outerS;
    stream.read(reinterpret_cast<char *>(&m), sizeof(int64_t));
    stream.read(reinterpret_cast<char *>(&n), sizeof(int64_t));
    stream.read(reinterpret_cast<char *>(&nnzS), sizeof(int64_t));
    stream.read(reinterpret_cast<char *>(&outerS), sizeof(int64_t));
    stream.read(reinterpret_cast<char *>(&innerS), sizeof(int64_t));

    this->referenceMask.resize(m, n);
    this->referenceMask.reserve(
        Eigen::VectorXi::Constant(
            this->variantGroups.size(), 
            this->sMax - this->sMin + 1
            )
        );
    this->referenceMask.makeCompressed();
    this->referenceMask.resizeNonZeros(nnzS);

    stream.read(reinterpret_cast<char *>(this->referenceMask.valuePtr()), sizeof(float) * nnzS);
    stream.read(reinterpret_cast<char *>(this->referenceMask.outerIndexPtr()), sizeof(int64_t) * outerS);
    stream.read(reinterpret_cast<char *>(this->referenceMask.innerIndexPtr()), sizeof(int64_t) * nnzS);
    this->referenceMask.finalize();

    // get variant allele names
    int nA;
    stream.read(reinterpret_cast<char *>(&nA), sizeof(int));
    std::unordered_map<std::string, int>().swap(this->variantAlleleNames);

    for (int i = 0; i < nA; ++i)
    {
        stream.read(reinterpret_cast<char *>(&l), sizeof(int));
        tempString = new char[l];
        stream.read(tempString, l);
        if (tempString[l-1] != '\0')
            throw std::runtime_error("");
        std::string a(tempString);
        this->variantAlleleNames[a] = i;
        delete[] tempString;
    }

    this->variantMask = std::vector<std::vector<Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t>>>(
        nA,
        std::vector<Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t>>(
            this->sMax - this->sMin + 1,
            Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t>()
        )
    );

    for (int i = 0; i < nA; ++i)
    {
        for (int s = this->sMin; s <= this->sMax; ++s)
        {
            stream.read(reinterpret_cast<char *>(&m), sizeof(int64_t));
            stream.read(reinterpret_cast<char *>(&n), sizeof(int64_t));
            stream.read(reinterpret_cast<char *>(&nnzS), sizeof(int64_t));
            stream.read(reinterpret_cast<char *>(&outerS), sizeof(int64_t));
            stream.read(reinterpret_cast<char *>(&innerS), sizeof(int64_t));

            this->variantMask[i][s - this->sMin].resize(m, n);
            this->variantMask[i][s - this->sMin].reserve(
                Eigen::VectorXi::Constant(this->variantGroups.size(), (this->variant.getAllJunctions().size() + 1) * (this->sMax - this->sMin + 1))
            );
            this->variantMask[i][s - this->sMin].makeCompressed();
            this->variantMask[i][s - this->sMin].resizeNonZeros(nnzS);

            stream.read(reinterpret_cast<char *>(this->variantMask[i][s - this->sMin].valuePtr()), sizeof(float) * nnzS);
            stream.read(reinterpret_cast<char *>(this->variantMask[i][s - this->sMin].outerIndexPtr()), sizeof(int64_t) * outerS);
            stream.read(reinterpret_cast<char *>(this->variantMask[i][s - this->sMin].innerIndexPtr()), sizeof(int64_t) * nnzS);
            this->variantMask[i][s - this->sMin].finalize();
        }
    }

    // free allocated char arrays
    stream.close();
}

int VariantProfile::getMinInsert()
{
    return this->sMin;
}

int VariantProfile::getMaxInsert()
{
    return this->sMax;
}

int VariantProfile::getMargin()
{
    return this->filterMargin;
}

int VariantProfile::getReadLength()
{
    return this->readLength;
}

bool VariantProfile::variantStructureIsPresent()
{
    return this->variantPresent;
}

bool VariantProfile::loadVariantStructure(std::string filename, std::string variantName)
{
    // load with variant parser
    try {
        variantParser parser(filename);
        int idx = -1;
        for (uint32_t i = 0; i < parser.getVariantNames().size(); ++i)
        {
            if (parser.getVariantNames()[i] == variantName)
            {
                idx = i;
                break;
            }
        }
        if (idx < 0)
        {
            std::cerr << "Could not find variant " << variantName << " in specified file (" << filename << ")." << std::endl;
            return false;
        }

        this->variant = complexVariant(parser.getVariantNames()[idx], parser.getAlleleNames()[idx], parser.getVariantJunctions()[idx], filename);
        createVariantChromosomeStructures();
    } 
    catch (std::runtime_error const& err)
    {
        std::cerr << "Could not read variant file (" << filename << ")" << std::endl;
        return false;
    }

    return true;
}

std::string VariantProfile::getName()
{
    return this->name;
}
