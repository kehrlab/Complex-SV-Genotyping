#include "allele.hpp"
#include "seqan/modifier/modifier_reverse.h"
#include "seqan/sequence/sequence_shortcuts.h"
#include "vMapManager.hpp"
#include <stdexcept>
#include <unordered_set>

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
    std::runtime_error(msg.c_str());
}

void Allele::clearSequences()
{
    this->chromosomeSequences = seqan::String<seqan::Dna5String>();
    this->sequenceIDs = seqan::String<seqan::CharString>();
}