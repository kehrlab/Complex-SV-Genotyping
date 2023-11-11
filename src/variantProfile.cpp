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
    const std::unordered_map<std::string, int> & contigInfos
) : variant(variant), filterMargin(filterMargin), overlap(overlap), readLength(readLength), sMin(sMin), sMax(sMax)
{
    this->variant.createAlleleMaps(this->filterMargin, contigInfos);
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

void VariantProfile::initMasks()
{
    initReferenceMask();
    initVariantMask();
}


std::vector<std::unordered_set<int>> VariantProfile::createIndexCombinations(std::vector<int> & firstIndices, std::vector<int> & secondIndices)
{
    std::vector<std::unordered_set<int>> indexGroups;

    std::unordered_set<int> currentGroup;
    for (int i = 0; i < firstIndices.size(); ++i)
    {
        currentGroup.insert(firstIndices[i]);
        indexGroups.push_back(currentGroup);

        std::unordered_set<int> tempGroup {currentGroup};
        for (int j = 1; j <= secondIndices.size(); ++j)
        {
            tempGroup.insert(secondIndices[secondIndices.size() - j]);
            indexGroups.push_back(tempGroup);

            std::unordered_set<int> tempGroup2 {currentGroup};
            for (int k = j; k <= secondIndices.size(); ++k)
            {
                tempGroup2.insert(secondIndices[secondIndices.size() - k]);
                indexGroups.push_back(tempGroup2);
            }
        }
    }
    return indexGroups;
}

std::vector<std::unordered_set<int>> VariantProfile::createIndexCombinations(std::vector<int> & firstIndices, std::vector<int> & secondIndices, std::unordered_map<int, std::unordered_set<int>> & junctionMatches)
{
    std::vector<std::unordered_set<int>> indexGroups;

    std::unordered_set<int> currentGroup;
    for (int i = 0; i < firstIndices.size(); ++i)
    {
        currentGroup.insert(firstIndices[i]);
        indexGroups.push_back(currentGroup);
        
        std::unordered_set<int> tempGroup {currentGroup};
        for (int j = 1; j <= secondIndices.size(); ++j)
        {
            tempGroup.insert(secondIndices[secondIndices.size() - j]);
            indexGroups.push_back(tempGroup);

            std::unordered_set<int> tempGroup2 {currentGroup};
            for (int k = j; k <= secondIndices.size(); ++k)
            {
                tempGroup2.insert(secondIndices[secondIndices.size() - k]);
                indexGroups.push_back(tempGroup2);
            }
        }
    }

    // create possible alternatives
    std::vector<std::unordered_set<int>> altIndexGroups;
    for (int i = 0; i < indexGroups.size(); ++i)
    {
        std::vector<std::unordered_set<int>> altGroups;
        std::vector<std::unordered_set<int>> tempGroups;

        std::vector<int> tempIndices;
        for (auto & idx : indexGroups[i])
            tempIndices.push_back(idx);
        
        altGroups.push_back(std::unordered_set<int>{tempIndices[0]});
        if (junctionMatches.find(tempIndices[0]) != junctionMatches.end())
        {
            for (auto & idx : junctionMatches[tempIndices[0]])
                altGroups.push_back(std::unordered_set<int>{idx});
        }

        for (int j = 1; j < tempIndices.size(); ++j)
        {
            tempGroups = altGroups;
            altGroups = std::vector<std::unordered_set<int>>();

            for (int k = 0; k < tempGroups.size(); ++k)
            {
                std::unordered_set<int> temp{tempGroups[k]};
                temp.insert(tempIndices[j]);
                altGroups.push_back(temp);

                if (junctionMatches.find(tempIndices[j]) != junctionMatches.end())
                {
                    for (auto & idx : junctionMatches[tempIndices[j]])
                    {
                        temp = std::unordered_set<int>{tempGroups[k]};
                        temp.insert(idx);
                        altGroups.push_back(temp);
                    }
                }
            }
        }
        for (auto & aG : altGroups)
            altIndexGroups.push_back(aG);
    }
    return altIndexGroups;
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
        std::vector<int> indices;
        std::vector<int> localIndices;
        std::vector<char> type;

        for (int i = 0; i < jRegion.breakpointIndices.size(); ++i) 
        {
            indices.push_back(jRegion.breakpointIndices[i]);
            localIndices.push_back(i);
            type.push_back('B');
        }
        for (int i = 0; i < jRegion.junctionIndices.size(); ++i)
        {
            indices.push_back(jRegion.junctionIndices[i]);
            localIndices.push_back(i);
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
        for (int i = 1; i < indices.size(); ++i)
        {
            int d = indices[i] - indices[i-1];
            if (type[i-1] == 'J')
                --d;

            if (d > (2*this->filterMargin))
            {
                ++regionIdx;
                grouping.push_back(std::vector<int>{i});
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
            int minIdx {300000000}, maxIdx {0};
            char minType {' '}, maxType{' '};
            for (auto & p : grouping[g])
            {
                int idx = indices[p];

                char & t = type[helper[p]];
                int & lI = localIndices[helper[p]];

                if (idx < minIdx)
                {
                    minType = t;
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
            int nextDistance{-1}, lastDistance{-1};
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
        int regionLength = 0;
        for (int i = 0; i < breakpoints.size(); ++i)
            if (breakpoints[i].getReferenceName() == chr)
                bpIndices.push_back(i);
        for (int i = 0; i < junctions.size(); ++i)
            if (junctions[i].getVariantRefName() == chr)
                njIndices.push_back(i);
        std::sort(bpIndices.begin(), bpIndices.end(), [&](const int & i, const int & j)
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
                    (isReverse ? 300000000 : junctions[njIndices[0]].getPositionLeft()),
                    isReverse
                }
            );
            regionLength = chrStructure.regions[0].getRegionEnd() - chrStructure.regions[0].getRegionStart();
            chrStructure.junctionIndices.push_back(regionLength);
            chrStructure.length = regionLength + 1;
            chrStructure.junctions.push_back(junctions[njIndices[0]]);
            // middle
            for (int j = 0; j < njIndices.size() - 1; ++j)
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
            chrStructure.regions.push_back(
                GenomicRegion {
                    junctions[tempIdx].getRefNameRight(),
                    (isReverse ? 0 : junctions[tempIdx].getPositionRight()),
                    (isReverse ? junctions[tempIdx].getPositionRight() : 300000000),
                    isReverse
                }
            );
            regionLength = chrStructure.regions[tempIdx].getRegionEnd() - chrStructure.regions[tempIdx].getRegionStart() + 1;
            chrStructure.length += regionLength;
        } 
        else 
        {
            chrStructure = JunctionRegion {
                std::vector<GenomicRegion> {
                    GenomicRegion{
                        chr,
                        0,
                        300000000,
                        false
                    }
                },
                std::vector<int>(),
                std::vector<int>(),
                std::vector<Junction>(),
                std::vector<Breakpoint>(),
                300000000
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
    for (int j = 0; j < jRegion.regions.size(); ++j)
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
            ++counter;
        }
    }
}

void VariantProfile::findPairAttributes(std::unordered_set<std::string> & groups, VariantRegions & vRegions)
{
    groups.insert("RF");
    for (int regionIdx = 0; regionIdx < vRegions.regions.size(); ++regionIdx)
    {
        auto & jRegion = vRegions.regions[regionIdx];
        for (int j = 0; j < jRegion.junctionIndices.size(); ++j)
        {
            // at every junction, try the smallest and largest possible insert size
            // and the nearest and furthest possible location
            // a pair may span more than one junction
            for (auto s : std::vector<int>{std::max(this->sMin, 2), this->sMax})
            {
                int beginIndexLeft = jRegion.junctionIndices[j] - (s - 2);
                int endIndexLeft = jRegion.junctionIndices[j] + 1;
                int endPosLeft = jRegion.junctions[j].getPositionRight();

                int beginIndexRight = jRegion.junctionIndices[j];
                int endIndexRight = jRegion.junctionIndices[j] + s - 1;
                int beginPosRight = jRegion.junctions[j].getPositionLeft();

                // consider junctions to the left
                for (int k = j; k >= 0; --k)
                {
                    int remaining = s - (jRegion.junctionIndices[j] - jRegion.junctionIndices[k] + 2);
                    bool fr {true};

                    if (remaining >= 0)
                    {   
                        // chromosomes
                        std::unordered_set<std::string> chromosomes;
                        chromosomes.insert(jRegion.junctions[k].getRefNameLeft());
                        chromosomes.insert(jRegion.junctions[j].getRefNameRight());
                        if (chromosomes.size() > 1)
                        {
                            std::string chromosomeString;
                            createIndexString(chromosomeString, chromosomes);
                            groups.insert(chromosomeString);
                        }

                        // orientation
                        if (jRegion.junctions[k].getDirectionLeft() < 0 && jRegion.junctions[j].getDirectionRight() < 0 && chromosomes.size() == 1)
                        {
                            fr = false;
                            groups.insert("FF");
                        }
                        else if (jRegion.junctions[k].getDirectionLeft() > 0 && jRegion.junctions[j].getDirectionRight() > 0 && chromosomes.size() == 1)
                        {
                            fr = false;
                            groups.insert("RR");
                        }

                        // the one clipped at the last region
                        int beginIdx = jRegion.junctionIndices[k] + 1;
                        int beginPos = jRegion.junctions[k].getPositionRight();
                        int sNew = endPosLeft - beginPos + 1;

                        if (chromosomes.size() == 1)
                        {
                            if (!fr)
                                sNew = std::abs(sNew);
                        } 
                        else
                        {
                            sNew = 0;
                        }
                        this->sMinMapped = std::min(this->sMinMapped, sNew);
                        this->sMaxMapped = std::max(this->sMaxMapped, sNew);

                        if (k > 0 && remaining >= (jRegion.junctionIndices[k] - jRegion.junctionIndices[k-1]))
                            continue;
                        if (k==0 && vRegions.distanceFromLast[regionIdx] > 0 && remaining >= (vRegions.distanceFromLast[regionIdx] + 2*this->filterMargin))
                            break;
                        
                        // the rest
                        beginIdx = beginIndexLeft;
                        beginPos = jRegion.junctions[k].getPositionLeft() +
                            jRegion.junctions[k].getDirectionLeft() *
                            remaining;
                        sNew = endPosLeft - beginPos + 1;

                        if (chromosomes.size() == 1)
                        {
                            if (!fr)
                                sNew = std::abs(sNew);
                        } 
                        else
                        {
                            sNew = 0;
                        }
                        this->sMinMapped = std::min(this->sMinMapped, sNew);
                        this->sMaxMapped = std::max(this->sMaxMapped, sNew);
                    }
                }

                // also on previous region(s)
                int tempIdx = regionIdx;
                int dTemp = jRegion.junctionIndices[j];
                while (vRegions.distanceFromLast[tempIdx] > 0)
                {
                    bool someLeft = true;
                    JunctionRegion & lastRegion = vRegions.regions[tempIdx - 1];
                    for (int k = lastRegion.junctionIndices.size() - 1; k >= 0; --k)
                    {
                        int d = vRegions.distanceFromLast[tempIdx] + dTemp + 
                            (lastRegion.length - 1 - lastRegion.junctionIndices[k]);
                        int remaining = s - (d+2);

                        bool fr {true};

                        if (remaining >= 0)
                        {   
                            // chromosomes
                            std::unordered_set<std::string> chromosomes;
                            chromosomes.insert(lastRegion.junctions[k].getRefNameLeft());
                            chromosomes.insert(jRegion.junctions[j].getRefNameRight());
                            if (chromosomes.size() > 1)
                            {
                                std::string chromosomeString;
                                createIndexString(chromosomeString, chromosomes);
                                groups.insert(chromosomeString);
                            }

                            // orientation
                            if (lastRegion.junctions[k].getDirectionLeft() < 0 && jRegion.junctions[j].getDirectionRight() < 0 && chromosomes.size() == 1)
                            {
                                fr = false;
                                groups.insert("FF");
                            }
                            else if (lastRegion.junctions[k].getDirectionLeft() > 0 && jRegion.junctions[j].getDirectionRight() > 0 && chromosomes.size() == 1)
                            {
                                fr = false;
                                groups.insert("RR");
                            }

                            // the one clipped at the last region
                            int beginIdx = lastRegion.junctionIndices[k] + 1;
                            int beginPos = lastRegion.junctions[k].getPositionRight();
                            int sNew = endPosLeft - beginPos + 1;

                            if (chromosomes.size() == 1)
                            {
                                if (!fr)
                                    sNew = std::abs(sNew);
                            } 
                            else
                            {
                                sNew = 0;
                            }
                            this->sMinMapped = std::min(this->sMinMapped, sNew);
                            this->sMaxMapped = std::max(this->sMaxMapped, sNew);

                            if (k > 0 && remaining >= (lastRegion.junctionIndices[k] - lastRegion.junctionIndices[k-1]))
                                continue;
                            
                            // the rest
                            beginPos = lastRegion.junctions[k].getPositionLeft() +
                                lastRegion.junctions[k].getDirectionLeft() *
                                remaining;
                            sNew = endPosLeft - beginPos + 1;

                            if (chromosomes.size() == 1)
                            {
                                if (!fr)
                                    sNew = std::abs(sNew);
                            } 
                            else
                            {
                                sNew = 0;
                            }
                            this->sMinMapped = std::min(this->sMinMapped, sNew);
                            this->sMaxMapped = std::max(this->sMaxMapped, sNew);
                        } else {
                            someLeft = false;
                            break;
                        }
                    }
                    if (!someLeft)
                        break;
                    
                    dTemp += (vRegions.distanceFromLast[tempIdx] + lastRegion.length - 1);
                    --tempIdx;
                }


                // consider junctions to the right
                for (int k = j; k < jRegion.junctionIndices.size(); ++k)
                {
                    int remaining = s - (jRegion.junctionIndices[k] - jRegion.junctionIndices[j] + 2);
                    bool fr {true};

                    if (remaining >= 0)
                    {   
                        // chromosomes
                        std::unordered_set<std::string> chromosomes;
                        chromosomes.insert(jRegion.junctions[j].getRefNameLeft());
                        chromosomes.insert(jRegion.junctions[k].getRefNameRight());
                        if (chromosomes.size() > 1)
                        {
                            std::string chromosomeString;
                            createIndexString(chromosomeString, chromosomes);
                            groups.insert(chromosomeString);
                        }

                        // orientation
                        if (jRegion.junctions[j].getDirectionLeft() < 0 && jRegion.junctions[k].getDirectionRight() < 0 && chromosomes.size() == 1)
                        {
                            fr = false;
                            groups.insert("FF");
                        }
                        else if (jRegion.junctions[j].getDirectionLeft() > 0 && jRegion.junctions[k].getDirectionRight() > 0 && chromosomes.size() == 1)
                        {
                            fr = false;
                            groups.insert("RR");
                        }

                        // the one clipped at the last region
                        int endIdx = jRegion.junctionIndices[k];
                        int endPos = jRegion.junctions[k].getPositionLeft();
                        int sNew = endPos - beginPosRight + 1;

                        if (chromosomes.size() == 1)
                        {
                            if (!fr)
                                sNew = std::abs(sNew);
                        } 
                        else
                        {
                            sNew = 0;
                        }
                        this->sMinMapped = std::min(this->sMinMapped, sNew);
                        this->sMaxMapped = std::max(this->sMaxMapped, sNew);

                        if (k < (jRegion.junctions.size() - 1) && remaining >= (jRegion.junctionIndices[k+1] - jRegion.junctionIndices[k]))
                            continue;
                        if (k==(jRegion.junctions.size() - 1) && vRegions.distanceToNext[regionIdx] > 0 && remaining >= (vRegions.distanceToNext[regionIdx] + 2*this->filterMargin))
                            break;
                        
                        // the rest
                        endIdx = endIndexRight;
                        endPos = jRegion.junctions[k].getPositionRight() +
                            jRegion.junctions[k].getDirectionRight() *
                            remaining;
                        sNew = endPos - beginPosRight + 1;

                        if (chromosomes.size() == 1)
                        {
                            if (!fr)
                                sNew = std::abs(sNew);
                        } 
                        else
                        {
                            sNew = 0;
                        }
                        this->sMinMapped = std::min(this->sMinMapped, sNew);
                        this->sMaxMapped = std::max(this->sMaxMapped, sNew);
                    }
                }

                // also on next region(s)
                tempIdx = regionIdx;
                dTemp = (jRegion.length - 1 - jRegion.junctionIndices[j]);
                while (vRegions.distanceToNext[tempIdx] > 0)
                {
                    bool someLeft = true;
                    JunctionRegion & nextRegion = vRegions.regions[tempIdx + 1];
                    for (int k = 0; k < nextRegion.junctionIndices.size(); ++k)
                    {
                        int d = vRegions.distanceToNext[tempIdx] + dTemp + 
                            nextRegion.junctionIndices[k];
                        int remaining = s - (d+2);

                        bool fr {true};

                        if (remaining >= 0)
                        {   
                            // chromosomes
                            std::unordered_set<std::string> chromosomes;
                            chromosomes.insert(nextRegion.junctions[k].getRefNameRight());
                            chromosomes.insert(jRegion.junctions[j].getRefNameLeft());
                            if (chromosomes.size() > 1)
                            {
                                std::string chromosomeString;
                                createIndexString(chromosomeString, chromosomes);
                                groups.insert(chromosomeString);
                            }

                            // orientation
                            if (nextRegion.junctions[k].getDirectionRight() < 0 && jRegion.junctions[j].getDirectionLeft() < 0 && chromosomes.size() == 1)
                            {
                                fr = false;
                                groups.insert("FF");
                            }
                            else if (nextRegion.junctions[k].getDirectionRight() > 0 && jRegion.junctions[j].getDirectionLeft() > 0 && chromosomes.size() == 1)
                            {
                                fr = false;
                                groups.insert("RR");
                            }

                            // the one clipped at the last region
                            int endIdx = nextRegion.junctionIndices[k];
                            int endPos = nextRegion.junctions[k].getPositionLeft();
                            int sNew = endPos - beginPosRight + 1;

                            if (chromosomes.size() == 1)
                            {
                                if (!fr)
                                    sNew = std::abs(sNew);
                            } 
                            else
                            {
                                sNew = 0;
                            }
                            this->sMinMapped = std::min(this->sMinMapped, sNew);
                            this->sMaxMapped = std::max(this->sMaxMapped, sNew);

                            if (k < nextRegion.junctions.size() - 1 && remaining >= (nextRegion.junctionIndices[k+1] - nextRegion.junctionIndices[k]))
                                continue;
                            
                            // the rest
                            endPos = nextRegion.junctions[k].getPositionRight() +
                                nextRegion.junctions[k].getDirectionRight() *
                                remaining;
                            sNew = endPos - beginPosRight + 1;

                            if (chromosomes.size() == 1)
                            {
                                if (!fr)
                                    sNew = std::abs(sNew);
                            } 
                            else
                            {
                                sNew = 0;
                            }
                            this->sMinMapped = std::min(this->sMinMapped, sNew);
                            this->sMaxMapped = std::max(this->sMaxMapped, sNew);
                        } else {
                            someLeft = false;
                            break;
                        }
                    }
                    if (!someLeft)
                        break;

                    dTemp += vRegions.distanceToNext[tempIdx] + nextRegion.length - 1;
                    ++tempIdx;
                }
            }
        }
    }
    // add a little buffer at the edges...
    this->sMinMapped -= 100;
    this->sMaxMapped += 100;
}

void VariantProfile::determineSpanningGroups(std::unordered_set<std::string> & groups, VariantRegions & vRegions)
{
    std::vector<std::unordered_set<int>> indexGroups;
    // consider all regions
    for (int v = 0; v < vRegions.regions.size(); ++v)
    {
        JunctionRegion & jRegion = vRegions.regions[v];
        // consider all breakpoints as possible starting positions for read pairs
        for (int i = 0; i < jRegion.breakpoints.size(); ++i)
        {
            std::vector<int> firstIndices;
            std::unordered_set<int> tempSet;

            // start first read as close to breakpoint as possible
            int firstReadEnd {jRegion.breakpointIndices[i] + this->readLength - this->overlap};
            firstIndices.push_back(jRegion.breakpoints[i].getID());
            tempSet.insert(jRegion.breakpoints[i].getID());
            indexGroups.push_back(tempSet);
            

            // check for other breakpoints that are spanned by the read
            for (int j = i + 1; j < jRegion.breakpoints.size(); ++j)
            {
                if (jRegion.breakpointIndices[j] <= firstReadEnd) 
                {
                    firstIndices.push_back(jRegion.breakpoints[j].getID());
                    tempSet.insert(jRegion.breakpoints[j].getID());
                    indexGroups.push_back(tempSet);
                }
                else
                {
                    break;
                }
            }

            // go over possible end positions for the second read pair
            int j = i + 1;
            for (; j < jRegion.breakpoints.size(); ++j)
            {
                if (jRegion.breakpointIndices[j] - jRegion.breakpointIndices[i] > this->sMax - 2*this->overlap)
                    break;
                
                std::vector<int> secondIndices;
                secondIndices.push_back(jRegion.breakpoints[j].getID());
                int secondReadStart {jRegion.breakpointIndices[j] - this->readLength + this->overlap};

                for (int k = j - 1; k > i; --k)
                {
                    if (jRegion.breakpointIndices[k] >= secondReadStart)
                        secondIndices.push_back(jRegion.breakpoints[k].getID());
                    else 
                        break;
                }

                std::vector<std::unordered_set<int>> tempGroups = createIndexCombinations(firstIndices, secondIndices);
                for (auto & set : tempGroups)
                    indexGroups.push_back(set);
            }

            // Sometimes it is necessary to consider breakpoints on next region(s) as well
            if (j != jRegion.breakpoints.size())
                continue;

            int tempIdx = v;
            int dTemp = jRegion.length - 1 - jRegion.breakpointIndices[i];
            while (vRegions.distanceToNext[tempIdx] > 0)
            {
                bool someLeft = true;
                JunctionRegion & nextRegion = vRegions.regions[tempIdx + 1];
                for (int j = 0; j < nextRegion.breakpoints.size(); ++j)
                {
                    int d = vRegions.distanceToNext[tempIdx] + dTemp + nextRegion.breakpointIndices[j];
                    if (d > this->sMax - 2*this->overlap)
                    {
                        someLeft = false;
                        break;
                    }

                    std::vector<int> secondIndices;
                    secondIndices.push_back(nextRegion.breakpoints[j].getID());
                    int secondReadStart {nextRegion.breakpointIndices[j] - this->readLength + this->overlap};

                    for (int k = j - 1; k >= 0; --k)
                    {
                        if (nextRegion.breakpointIndices[k] >= secondReadStart)
                            secondIndices.push_back(nextRegion.breakpoints[k].getID());
                        else 
                            break;
                    }

                    std::vector<std::unordered_set<int>> tempGroups = createIndexCombinations(firstIndices, secondIndices);
                    for (auto & set : tempGroups)
                        indexGroups.push_back(set);
                }
                if (!someLeft)
                    break;
                dTemp += vRegions.distanceToNext[tempIdx] + nextRegion.length - 1;
                ++tempIdx;
            }
        }
    }

    // create final unique group strings from breakpoint indices
    std::string spanningString;
    for (int i = 0; i < indexGroups.size(); ++i)
    {
        createIndexString(spanningString, indexGroups[i]);
        spanningString = "spanning_" + spanningString;
        groups.insert(spanningString);
    }
}

std::unordered_map<int, std::unordered_set<int>> VariantProfile::getJunctionAmbiguities(VariantRegions & vRegions)
{
    std::unordered_map<int, std::unordered_set<int>> junctionMatches;
    for (auto & jRegion : vRegions.regions)
    {
        for (int i = 0; i < jRegion.junctions.size(); ++i)
        {
            Junction & j = jRegion.junctions[i];
            std::unordered_set<int> indices;
            
            for (auto & tempRegion : vRegions.regions)
            {
                for (int k = 0; k < tempRegion.junctions.size(); ++k)
                {
                    if (j.getID() == tempRegion.junctions[k].getID())
                        continue;
                    Junction & j2 = tempRegion.junctions[k];
                    
                    if (j.getRefNameLeft() == j2.getRefNameLeft() && std::abs(j.getPositionLeft() - j2.getPositionLeft()) <= 20)
                        indices.insert(j2.getID());
                    else if (j.getRefNameLeft() == j2.getRefNameRight() && std::abs(j.getPositionLeft() - j2.getPositionRight()) <= 20)
                        indices.insert(j2.getID());
                    else if (j.getRefNameRight() == j2.getRefNameRight() && std::abs(j.getPositionRight() - j2.getPositionRight()) <= 20)
                        indices.insert(j2.getID());
                    else if (j.getRefNameRight() == j2.getRefNameLeft() && std::abs(j.getPositionRight() - j2.getPositionLeft()) <= 20)
                        indices.insert(j2.getID());
                }
            }
            if (indices.size() > 0)
                junctionMatches[j.getID()] = indices;
        }
    }

    return junctionMatches;
}

void VariantProfile::determineSplitGroups(std::unordered_set<std::string> & groups, VariantRegions & vRegions)
{
    std::vector<std::unordered_set<int>> indexGroups;
    std::unordered_map<int, std::unordered_set<int>> junctionMatches = getJunctionAmbiguities(vRegions);

    // consider all regions
    for (int v = 0; v < vRegions.regions.size(); ++v)
    {
        JunctionRegion & jRegion = vRegions.regions[v];
        // consider all junctions as possible starting positions for read pairs
        for (int i = 0; i < jRegion.junctions.size(); ++i)
        {
            std::vector<int> firstIndices;
            std::unordered_set<int> tempSet;

            // start first read as close to junction as possible
            int firstReadEnd {jRegion.junctionIndices[i] + this->readLength - this->overlap};
            firstIndices.push_back(jRegion.junctions[i].getID());
            tempSet.insert(jRegion.junctions[i].getID());
            indexGroups.push_back(tempSet);

            // check for other junctions that are spanned by the read
            for (int j = i + 1; j < jRegion.junctions.size(); ++j)
            {
                if (jRegion.junctionIndices[j] + 1 <= firstReadEnd)
                {
                    firstIndices.push_back(jRegion.junctions[j].getID());
                    tempSet.insert(jRegion.junctions[j].getID());
                    indexGroups.push_back(tempSet);
                }
                else
                {
                    break;
                }
            }

            // go over possible end positions for the second read in pair
            int j = i + 1;
            
            for (; j < jRegion.junctions.size(); ++j)
            {
                if (jRegion.junctionIndices[j] + 1 - jRegion.junctionIndices[i] > this->sMax - 2*this->overlap)
                    break;
                
                std::vector<int> secondIndices;
                secondIndices.push_back(jRegion.junctions[j].getID());
                int secondReadStart {jRegion.junctionIndices[j] - this->readLength + this->overlap};

                for (int k = j - 1; k > i; --k)
                {
                    if (jRegion.junctionIndices[k] >= secondReadStart)
                        secondIndices.push_back(jRegion.junctions[k].getID());
                    else
                        break;
                }
		
                std::vector<std::unordered_set<int>> tempGroups = createIndexCombinations(firstIndices, secondIndices, junctionMatches);
                for (auto & set : tempGroups)
                    indexGroups.push_back(set);
            }

            // Sometimes it is necessary to consider junctions on next region as well
            if (j != jRegion.junctions.size())
                continue;

            int tempIdx = v;
            int dTemp = jRegion.length - 1 - jRegion.junctionIndices[i];
            while (vRegions.distanceToNext[tempIdx] > 0)
            {
                bool someLeft = true;
                JunctionRegion & nextRegion = vRegions.regions[v+1];
                for (int j = 0; j < nextRegion.junctions.size(); ++j)
                {
                    int d = vRegions.distanceToNext[tempIdx] + dTemp + nextRegion.junctionIndices[j];
                    if (d > this->sMax - 2*this->overlap)
                    {
                        someLeft = false;
                        break;
                    }

                    std::vector<int> secondIndices;
                    secondIndices.push_back(nextRegion.junctions[j].getID());
                    int secondReadStart {nextRegion.junctionIndices[j] - this->readLength + this->overlap};

                    for (int k = j - 1; k >= 0; --k)
                    {
                        if (nextRegion.junctionIndices[k] >= secondReadStart) 
                            secondIndices.push_back(nextRegion.junctions[k].getID());
                        else
                            break;
                    }
		
                    std::vector<std::unordered_set<int>> tempGroups = createIndexCombinations(firstIndices, secondIndices, junctionMatches);
                    for (auto & set : tempGroups)
                        indexGroups.push_back(set);
                }

                dTemp += vRegions.distanceToNext[tempIdx] + nextRegion.length - 1;
                ++tempIdx;
            }
        }
    }

    // create final unique group strings from junction indices
    std::string splitString;
    for (int i = 1; i <= indexGroups.size(); ++i)
    {
        createIndexString(splitString, indexGroups[i-1]);
        splitString = "split_" + splitString;
        groups.insert(splitString);
    }
    groups.insert("split_ambiguous");
}

void VariantProfile::initReferenceMask()
{
    this->referenceMask = Eigen::SparseMatrix<float, Eigen::RowMajor>(this->sMaxMapped - this->sMinMapped + 1, this->variantGroups.size());
}

void VariantProfile::initVariantMask()
{
    for (int i = 0; i < this->variant.getAlleles().size() - 1; ++i)
    {
        this->variantMask.push_back(
            std::vector<Eigen::SparseMatrix<float, Eigen::RowMajor>>(
                this->sMax - this->sMin + 1, 
                Eigen::SparseMatrix<float, Eigen::RowMajor>(this->sMaxMapped - this->sMinMapped + 1, this->variantGroups.size())
            )
        );
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

            for (int j = 0; j < maps.size(); ++j)
            {
                VariantMap & map = maps[j];
                
                int mapLength = map.totalLength;
                
                for (unsigned i = 0; i < mapLength; ++i)
                {
                    for (int s = this->sMax; s >= this->sMin; --s)
                    {
                        if (s >= mapLength || i+s >= mapLength)
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
}

inline void VariantProfile::addSimulatedTemplateToMask(int & status, VariantMap & map, int i, int s, Allele & allele)
{
    ReadTemplate sT = map.simulateTemplate(i, s, this->readLength);
    if (!sT.isInterestingReadPair(this->filter))
    {
        status = -1;
        return;
    }

    sT.findSpanningReads(this->variant.getAllBreakpoints());
    if (sT.containsSuspectedSplit())
        sT.findSplitReads(allele.getNovelJunctions(), this->chromosomeStructures);
    sT.determineLocationStrings();

    std::string orientation = sT.getOrientation();
    int insertSize = sT.getInsertSize();
    bool split = sT.containsSplitRead();
    bool spanning = sT.containsSpanningRead();
    bool interChromosome = sT.alignsAcrossChromosomes();
    std::string junctionString = sT.getJunctionString();
    std::string breakpointString = sT.getBreakpointString();
    std::string chrString = sT.getChromosomeString();
    
    status = 1;
    addValueToMask(allele, s, insertSize, orientation, junctionString, breakpointString, chrString);
    return;
}

inline void VariantProfile::addValueToMask(Allele & allele, int sOld, int sNew, std::string & orientation, std::string & jString, std::string & bpString, std::string & chromosomes)
{
    std::string group = GenotypeDistribution::determineGroup(orientation, jString, bpString, chromosomes);
    
    if (group == "" || this->variantGroups.find(group) == this->variantGroups.end()) {
        std::cout << "Warning: Could not assign read pair to any group. May cause problems if too frequent." << std::endl;
        std::cout << "Variant: " << this->variant.getName() << std::endl; 
        std::cout << "Group: " << group << std::endl;
        return;
    }

    if (sNew < this->sMinMapped || sNew >= this->sMaxMapped)
	    std::cerr << "Insert size error. Profile boundaries: " << this->sMinMapped << "\t" << this->sMaxMapped << "\tMapped insert size: " << sNew << std::endl;

    int gIdx = this->variantGroups[group];

    if (allele.getName() == "REF")
        this->referenceMask.coeffRef(sOld-this->sMinMapped, gIdx) += 1;
    else
        this->variantMask[this->variantAlleleNames[allele.getName()]][sOld-this->sMin].coeffRef(sNew - this->sMinMapped, gIdx) += 1;
    return;
}

const Eigen::SparseMatrix<float, Eigen::RowMajor> & VariantProfile::getVariantMask(int s)
{
    return this->variantMask[0][s];
}

const Eigen::SparseMatrix<float, Eigen::RowMajor> & VariantProfile::getVariantMask(const std::string & alleleName, int s)
{
    return this->variantMask[this->variantAlleleNames[alleleName]][s];
}

const Eigen::SparseMatrix<float, Eigen::RowMajor> & VariantProfile::getReferenceMask()
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
    std::vector<Eigen::SparseMatrix<float, Eigen::RowMajor>> variantAlleleDists(
        alleleNames.size(), 
        Eigen::SparseMatrix<float, Eigen::RowMajor>(this->sMaxMapped - this->sMinMapped + 1, this->variantGroups.size())
    );

    for (int s = this->sMin; s <= this->sMax; ++s)
    {
        float p = libraryDistribution.getProbability(s);
        for (int v = 0; v < alleleNames.size(); ++v)
        {
            if (alleleNames[v] == "REF")
                variantAlleleDists[0].row(s - this->sMinMapped) = this->referenceMask.row(s-this->sMinMapped) * p;
            else
                variantAlleleDists[v] += this->variantMask[this->variantAlleleNames[alleleNames[v]]][s - this->sMin] * p;
        }
    }

    // mix distributions
    float majorFactor = 1. - eps;
    float minorFactor = eps / variantAlleleDists.size();

    std::vector<Eigen::SparseMatrix<float, Eigen::RowMajor>> variantGtDists;
    std::vector<std::string> gtNames;

    Eigen::SparseMatrix<float, Eigen::RowMajor> temp(variantAlleleDists[0].rows(), variantAlleleDists[0].cols());
    for (int i = 0; i < variantAlleleDists.size(); ++i)
    {
        temp = 2 * variantAlleleDists[i] * majorFactor;
        for (int j = 0; j < variantAlleleDists.size(); ++j)
        {
            if (j == i)
                continue;
            temp += minorFactor * variantAlleleDists[j];
        }
        variantGtDists.push_back(temp);
        gtNames.push_back(alleleNames[i] + "/" + alleleNames[i]);

        for (int j = i + 1; j < variantAlleleDists.size(); ++j)
        {
            temp  = (variantAlleleDists[i] + variantAlleleDists[j]) * majorFactor;
            for (int k = 0; k < variantAlleleDists.size(); ++k)
            {
                if (k == i || k == j)
                    continue;
                temp += minorFactor * variantAlleleDists[k];
            }
            variantGtDists.push_back(temp);
            gtNames.push_back(alleleNames[i] + "/" + alleleNames[j]);
        }
    }

    // scale distributions
    for (auto & dist : variantGtDists)
    {
        float sum = dist.sum();
        dist /= sum;
	// std::cout << "0 / chr20chr21: " << dist.coeffRef(0 - this->sMinMapped, this->variantGroups["chr20chr21"]) << std::endl;
    }
    

    // create GenotypeDistributions and map
    float minProb = 1;
    for (int i = 0; i < variantGtDists.size(); ++i) 
    {
        distributions[gtNames[i]] = GenotypeDistribution(variantGtDists[i], this->variantGroups, this->sMinMapped, this->sMaxMapped);
        if (distributions[gtNames[i]].getMinProbability() > 0)
            minProb = std::min(minProb, distributions[gtNames[i]].getMinProbability());
    }
    for (int i = 0; i < variantGtDists.size(); ++i)
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
    stream.write(reinterpret_cast<const char *>(&this->sMinMapped), sizeof(int));
    stream.write(reinterpret_cast<const char *>(&this->sMaxMapped), sizeof(int));

    //if (this->variant.getName() == "Translocation_0")
    ///	     std::cout << "Writing...\tsMin: " << this->sMinMapped << "\t" << this->sMaxMapped << std::endl;

    // write filter margin
    stream.write(reinterpret_cast<const char *>(&this->filterMargin), sizeof(int));

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
    this->referenceMask.makeCompressed()  ;
    int m      = referenceMask.rows()     ;
    int n      = referenceMask.cols()     ;
    int nnzS   = referenceMask.nonZeros() ;
    int outerS = referenceMask.outerSize();
    int innerS = referenceMask.innerSize();
    stream.write(reinterpret_cast<const char *>(&m), sizeof(int));
    stream.write(reinterpret_cast<const char *>(&n), sizeof(int));
    stream.write(reinterpret_cast<const char *>(&nnzS), sizeof(int));
    stream.write(reinterpret_cast<const char *>(&outerS), sizeof(int));
    stream.write(reinterpret_cast<const char *>(&innerS), sizeof(int));

    stream.write(reinterpret_cast<const char *>(referenceMask.valuePtr()), sizeof(float)*nnzS);
    stream.write(reinterpret_cast<const char *>(referenceMask.outerIndexPtr()), sizeof(int)*outerS);
    stream.write(reinterpret_cast<const char *>(referenceMask.innerIndexPtr()), sizeof(int)*nnzS);

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
    for (int i = 0; i < this->variantMask.size(); ++i)
    {
        for (int s = 0; s < this->variantMask[i].size(); ++s)
        {
            this->variantMask[i][s].makeCompressed();

            int m      = this->variantMask[i][s].rows()     ;
            int n      = this->variantMask[i][s].cols()     ;
            int nnzS   = this->variantMask[i][s].nonZeros() ;
            int outerS = this->variantMask[i][s].outerSize();
            int innerS = this->variantMask[i][s].innerSize();

            stream.write(reinterpret_cast<const char *>(&m), sizeof(int))     ;
            stream.write(reinterpret_cast<const char *>(&n), sizeof(int))     ;
            stream.write(reinterpret_cast<const char *>(&nnzS), sizeof(int))  ;
            stream.write(reinterpret_cast<const char *>(&outerS), sizeof(int));
            stream.write(reinterpret_cast<const char *>(&innerS), sizeof(int));

            stream.write(reinterpret_cast<const char *>(this->variantMask[i][s].valuePtr()), sizeof(float)*nnzS)     ;
            stream.write(reinterpret_cast<const char *>(this->variantMask[i][s].outerIndexPtr()), sizeof(int)*outerS);
            stream.write(reinterpret_cast<const char *>(this->variantMask[i][s].innerIndexPtr()), sizeof(int)*nnzS)  ;
        }
    }
    // std::cout << "m, n: " << this->variantMask[0][0].rows() << ", " << this->variantMask[0][0].cols() << std::endl;
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

    // try load the actual variant
    this->variantPresent = loadVariantStructure(variantPath, this->name);

    // get read length, insert minima and maxima and filter margin
    stream.read(reinterpret_cast<char *>(&this->readLength), sizeof(int));
    stream.read(reinterpret_cast<char *>(&this->sMin), sizeof(int));
    stream.read(reinterpret_cast<char *>(&this->sMax), sizeof(int));
    stream.read(reinterpret_cast<char *>(&this->sMinMapped), sizeof(int));
    stream.read(reinterpret_cast<char *>(&this->sMaxMapped), sizeof(int));
    stream.read(reinterpret_cast<char *>(&this->filterMargin), sizeof(int));

    if (variantPresent)
    {
        this->variant.setFilterMargin(this->filterMargin);
        this->filter = ReadPairFilter(this->variant.getAllBreakpoints(), this->filterMargin, 0);
    }
    //if (this->variant.getName() == "Translocation_0")
    //    std::cout << "Reading...\tsMin: " << this->sMinMapped << "\t" << this->sMaxMapped << std::endl;

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
    int m, n, nnzS, innerS, outerS;
    stream.read(reinterpret_cast<char *>(&m), sizeof(int));
    stream.read(reinterpret_cast<char *>(&n), sizeof(int));
    stream.read(reinterpret_cast<char *>(&nnzS), sizeof(int));
    stream.read(reinterpret_cast<char *>(&outerS), sizeof(int));
    stream.read(reinterpret_cast<char *>(&innerS), sizeof(int));

    this->referenceMask.resize(m, n);
    this->referenceMask.makeCompressed();
    this->referenceMask.resizeNonZeros(nnzS);

    stream.read(reinterpret_cast<char *>(this->referenceMask.valuePtr()), sizeof(float) * nnzS);
    stream.read(reinterpret_cast<char *>(this->referenceMask.outerIndexPtr()), sizeof(int) * outerS);
    stream.read(reinterpret_cast<char *>(this->referenceMask.innerIndexPtr()), sizeof(int) * nnzS);
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

    // read the variant allele masks
    this->variantMask = std::vector<std::vector<Eigen::SparseMatrix<float, Eigen::RowMajor>>>(
        nA,
        std::vector<Eigen::SparseMatrix<float, Eigen::RowMajor>>(
            this->sMax - this->sMin + 1,
            Eigen::SparseMatrix<float, Eigen::RowMajor>(this->sMaxMapped - this->sMinMapped + 1, this->variantGroups.size())
        )
    );

    for (int i = 0; i < nA; ++i)
    {
        for (int s = this->sMin; s <= this->sMax; ++s)
        {
            stream.read(reinterpret_cast<char *>(&m), sizeof(int));
            stream.read(reinterpret_cast<char *>(&n), sizeof(int));
            stream.read(reinterpret_cast<char *>(&nnzS), sizeof(int));
            stream.read(reinterpret_cast<char *>(&outerS), sizeof(int));
            stream.read(reinterpret_cast<char *>(&innerS), sizeof(int));

            this->variantMask[i][s - this->sMin].resize(m, n);
            this->variantMask[i][s - this->sMin].makeCompressed();
            this->variantMask[i][s - this->sMin].resizeNonZeros(nnzS);

            stream.read(reinterpret_cast<char *>(this->variantMask[i][s - this->sMin].valuePtr()), sizeof(float) * nnzS);
            stream.read(reinterpret_cast<char *>(this->variantMask[i][s - this->sMin].outerIndexPtr()), sizeof(int) * outerS);
            stream.read(reinterpret_cast<char *>(this->variantMask[i][s - this->sMin].innerIndexPtr()), sizeof(int) * nnzS);
            this->variantMask[i][s - this->sMin].finalize();
        }
    }
    // std::cout << "m, n: " << this->variantMask[0][0].rows() << ", " << this->variantMask[0][0].cols() << std::endl;

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
        for (int i = 0; i < parser.getVariantNames().size(); ++i)
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
    } 
    catch (std::runtime_error)
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
