#include "allele.hpp"

Allele::Allele()
{}

Allele::Allele(std::vector<Junction> novelJunctions)
{
    this->alleleName = "";
    this->novelJunctions = novelJunctions;
    convertJunctionsToBreakpoints();
}

Allele::Allele(std::vector<Breakpoint> breakpoints)
{
    this->alleleName = "REF";
    this->referenceBreakpoints = breakpoints;
    removeDuplicateBreakpoints();
    reorderBreakPoints();
    joinBreakpoints();
}

Allele::Allele(std::vector<Junction> novelJunctions, std::string alleleName)
{
    this->alleleName = alleleName;
    this->novelJunctions = novelJunctions;
    convertJunctionsToBreakpoints();
}

void Allele::convertJunctionsToBreakpoints()
{
    splitAllJunctionsToBreakpoints();
    removeDuplicateBreakpoints();
    reorderBreakPoints();
    joinBreakpoints();
}

void Allele::splitAllJunctionsToBreakpoints()
{
    unsigned id = 1;
    for (auto junction : this->novelJunctions)
    {
        if (junction.hasValidLeftSide())
        {
            Breakpoint leftBp = junction.leftSideToBreakpoint(id);
            this->referenceBreakpoints.push_back(leftBp);
            ++id;
        }
        if (junction.hasValidRightSide())
        {
            Breakpoint rightBp = junction.rightSideToBreakpoint(id);
            this->referenceBreakpoints.push_back(rightBp);
            ++id;
        }  
    }
}

void Allele::removeDuplicateBreakpoints()
{
    for (auto it = this->referenceBreakpoints.begin(); it != this->referenceBreakpoints.end(); ++it)
    {
        for (auto it1 = it + 1; it1 != this->referenceBreakpoints.end(); )
        {
            if (*it1 == *it)
                it1 = this->referenceBreakpoints.erase(it1);
            else
                ++it1;
        }
    }
}

void Allele::reorderBreakPoints()
{
    std::unordered_set<std::string> refNames = getUniqueBpRefNames();
    for (auto refName : refNames)
        sortBreakpointsOnReference(refName);
}

void Allele::joinBreakpoints()
{
	for (auto bp = this->referenceBreakpoints.begin(); bp != this->referenceBreakpoints.end();++bp)
	{
		for (auto bp1 = bp + 1; bp1 != this->referenceBreakpoints.end(); ++bp1)
        {
            if (bp->getReferenceName() == bp1->getReferenceName() && std::abs(bp->getPosition() - bp1->getPosition()) < 2)
            {
                this->referenceBreakpoints.erase(bp1);
                break;
            }
        }
	}
}

std::unordered_set<std::string> Allele::getUniqueBpRefNames()
{
    std::unordered_set<std::string> refNames;
    for (auto bp : this->referenceBreakpoints)
        refNames.insert(bp.getReferenceName());
    return refNames;
}

void Allele::sortBreakpointsOnReference(std::string refName)
{
    std::vector<Breakpoint>::iterator firstElement, lastElement;

    for (auto it = this->referenceBreakpoints.begin(); it != this->referenceBreakpoints.end();)
    {
        if (it->getReferenceName() == refName)
        {
            firstElement = it;
            while (it->getReferenceName() == refName)
            {
                lastElement = ++it;
                if (it == this->referenceBreakpoints.end())
                    break;
            }
            break;
        }
        ++it;
    }

    std::sort(firstElement, lastElement, Breakpoint::compareBreakpoints);
}

std::vector<Junction> & Allele::getNovelJunctions()
{
    return this->novelJunctions;
}

std::vector<Breakpoint> & Allele::getReferenceBreakpoints()
{
    return this->referenceBreakpoints;
}

std::string Allele::getName()
{
    return this->alleleName;
}

void Allele::print()
{
    std::cout << this->alleleName << ":" << std::endl;
    std::cout << "----------------------" << std::endl;
    std::cout << "REFERENCE BREAKPOINTS:" << std::endl;
    std::cout << "----------------------" << std::endl;
    for (auto it : this->referenceBreakpoints)
        it.print();
    std::cout << "----------------------" << std::endl;
    std::cout << "   NOVEL JUNCTIONS:   " << std::endl;
    std::cout << "----------------------" << std::endl;
    for (auto it : this->novelJunctions)
        it.print();
    std::cout << "----------------------" << std::endl;
    std::cout << "   Chromosome Maps:   " << std::endl;
    std::cout << "----------------------" << std::endl;
    for (unsigned i = 0; i < this->chromosomeMaps.size(); ++i)
    {
        std::cout << this->chromosomeNames[i] << std::endl;
        this->chromosomeMaps[i].print();
    }
    std::cout << "----------------------" << std::endl;
}

void Allele::createChromosomeMaps(std::vector<Breakpoint> & allBreakpoints, int filterMargin, std::unordered_map<std::string, int> contigLengths)
{
    std::unordered_set<std::string> cNames;
    for (Breakpoint bp : allBreakpoints)
        cNames.insert(bp.getReferenceName());
    for (std::string cName : cNames) 
        createChromosomeMap(cName, allBreakpoints, filterMargin, contigLengths);
}

void Allele::createChromosomeMap(std::string cName, std::vector<Breakpoint> & allBreakpoints, int filterMargin, std::unordered_map<std::string, int> contigLengths)
{
    std::vector<Junction> chromosomeJunctions = getJunctionsOnChromosome(cName);
    VariantMapManager chromosomeMap(cName, filterMargin, contigLengths);
    if (chromosomeJunctions.size() > 0)
        chromosomeMap.createMapFromJunctions(chromosomeJunctions);
    chromosomeMap.createBreakpointMaps(allBreakpoints);

    this->chromosomeMaps.push_back(chromosomeMap);
    this->chromosomeNames.push_back(cName);
}

std::vector<Junction> Allele::getJunctionsOnChromosome(std::string cName)
{
    std::vector<Junction> junctions;
    for (auto junction : getNovelJunctions())
        if (junction.getVariantRefName() == cName)
            junctions.push_back(junction);
    return junctions;
}

void Allele::createAlleleSequences(SeqFileHandler & seqFileHandler)
{
    if (!seqFileHandler.isOpen())
        return;
    
    for (unsigned i = 0; i < this->chromosomeMaps.size(); ++i)
    {
        seqan::String<seqan::Dna5String> chrSequences = chromosomeMaps[i].getChromosomeSequences(seqFileHandler);
        for (unsigned j = 0; j < seqan::length(chrSequences); ++j)
        {
            seqan::CharString id = (chromosomeNames[i] + "_" + std::to_string(j)).c_str();
            seqan::appendValue(this->chromosomeSequences, chrSequences[j]);
            seqan::appendValue(this->sequenceIDs, id);
        }
    }
}

seqan::String<seqan::CharString> Allele::getSequenceIDs()
{
    return this->sequenceIDs;
}

seqan::String<seqan::Dna5String> Allele::getChromosomeSequences()
{
    return this->chromosomeSequences;
}

void Allele::writeSequence(SeqFileHandler & seqFileHandler, std::string directory)
{
    std::string filename = directory + "/" + this->alleleName + ".fa";
    seqFileHandler.writeSequencesToFile(this->sequenceIDs, this->chromosomeSequences, filename);
}

std::vector<VariantMapManager> & Allele::getChromosomeMaps()
{
    return this->chromosomeMaps;
}

std::vector<std::string> & Allele::getChromosomeNames()
{
    return this->chromosomeNames;
}

VariantMapManager & Allele::getChromosomeMap(std::string cName)
{
    for (unsigned i = 0; i < this->chromosomeNames.size(); ++i)
        if (this->chromosomeNames[i] == cName)
            return this->chromosomeMaps[i];

    std::string msg = "Could not find map for chromosome " + cName;
    throw std::runtime_error(msg.c_str());
}

void Allele::clearSequences()
{
    this->chromosomeSequences = seqan::String<seqan::Dna5String>();
    this->sequenceIDs = seqan::String<seqan::CharString>();
}


std::unordered_map<std::string, JunctionRegion> Allele::getChromosomeStructures(std::vector<Breakpoint> & breakpoints)
{
    std::unordered_map<std::string, JunctionRegion> structures;
    // find unique chromosomes
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
        for (int i = 0; i < this->novelJunctions.size(); ++i)
            if (this->novelJunctions[i].getVariantRefName() == chr)
                njIndices.push_back(i);
        std::sort(bpIndices.begin(), bpIndices.end(), [&](const int & i, const int & j)
        {
            return breakpoints[i].getPosition() < breakpoints[j].getPosition();
        });

        if (njIndices.size() > 0)
        {
            // create a set of regions describing the chromosome structure
            // left
            bool isReverse {this->novelJunctions[njIndices[0]].getDirectionLeft() > 0};
            chrStructure.regions.push_back(
                GenomicRegion {
                    this->novelJunctions[njIndices[0]].getRefNameLeft(),
                    (isReverse ? this->novelJunctions[njIndices[0]].getPositionLeft() : 0),
                    (isReverse ? 300000000 : this->novelJunctions[njIndices[0]].getPositionLeft()),
                    isReverse
                }
            );
            regionLength = chrStructure.regions[0].getRegionEnd() - chrStructure.regions[0].getRegionStart();
            chrStructure.junctionIndices.push_back(regionLength);
            chrStructure.length = regionLength + 1;
            chrStructure.junctions.push_back(this->novelJunctions[njIndices[0]]);
            // middle
            for (int j = 0; j < njIndices.size() - 1; ++j)
            {
                if (this->novelJunctions[njIndices[j]].getRefNameRight() != this->novelJunctions[njIndices[j + 1]].getRefNameLeft())
                    throw std::runtime_error("ERROR: RefNames of adjacent junctions do not match!");
                if (this->novelJunctions[njIndices[j]].getDirectionRight() == this->novelJunctions[njIndices[j + 1]].getDirectionLeft())
                    throw std::runtime_error("ERROR: Directions of adjacent junctions are incompatible!");
                isReverse = (this->novelJunctions[njIndices[j]].getDirectionRight() < 0);
                chrStructure.regions.push_back(
                    GenomicRegion {
                        this->novelJunctions[njIndices[j]].getRefNameRight(),
                        (isReverse ? this->novelJunctions[njIndices[j+1]].getPositionLeft() : this->novelJunctions[njIndices[j]].getPositionRight()),
                        (isReverse ? this->novelJunctions[njIndices[j]].getPositionRight() : this->novelJunctions[njIndices[j+1]].getPositionLeft()),
                        isReverse
                    }
                );
                regionLength = chrStructure.regions[j+1].getRegionEnd() - chrStructure.regions[j+1].getRegionStart() + 1;
                chrStructure.junctionIndices.push_back(
                    regionLength + chrStructure.junctionIndices[j]
                );
                chrStructure.length += regionLength;
                chrStructure.junctions.push_back(this->novelJunctions[njIndices[j+1]]);
            }
            // right
            int tempIdx {njIndices[njIndices.size() - 1]};
            isReverse = this->novelJunctions[tempIdx].getDirectionRight() < 0;
            chrStructure.regions.push_back(
                GenomicRegion {
                    this->novelJunctions[tempIdx].getRefNameRight(),
                    (isReverse ? 0 : this->novelJunctions[tempIdx].getPositionRight()),
                    (isReverse ? this->novelJunctions[tempIdx].getPositionRight() : 300000000),
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
        
        structures[chr] = chrStructure;
    }
    return structures;
}

void Allele::insertBreakpoint(Breakpoint & bp, JunctionRegion & jRegion)
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