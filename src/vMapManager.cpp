#include "vMapManager.hpp"

VariantMapManager::VariantMapManager()
{}

VariantMapManager::VariantMapManager(std::string chromosome, int filterMargin, std::unordered_map<std::string, int> contigLengths)
{
    this->variantMapPresent = false;
    setChromosomeName(chromosome);
    setFilterMargin(filterMargin);
    setContigLengths(contigLengths);
}

VariantMapManager::VariantMapManager(std::string chromosome, std::vector<Junction> novelJunctions, std::vector<Breakpoint> breakpoints, int filterMargin, std::unordered_map<std::string, int> contigLengths)
{
    this->variantMapPresent = false;
    setChromosomeName(chromosome);
    setFilterMargin(filterMargin);
    setContigLengths(contigLengths);
    createMapFromJunctions(novelJunctions);
    createBreakpointMaps(breakpoints);
}

void VariantMapManager::setChromosomeName(std::string chromosome)
{
    this->chromosome = chromosome;
}

void VariantMapManager::setFilterMargin(int filterMargin)
{
    this->filterMargin = filterMargin;
}

void VariantMapManager::setContigLengths(std::unordered_map<std::string, int> contigLengths)
{
    this->contigLengths = contigLengths;
}

void VariantMapManager::createMapFromJunctions(std::vector<Junction> novelJunctions)
{
    VariantMap map;
    fillMap(map, novelJunctions);
    this->variantMap = map;
    this->variantMapPresent = true;
    createReferenceRegions();
}

void VariantMapManager::createBreakpointMaps(std::vector<Breakpoint> breakpoints)
{
    for (Breakpoint bp : breakpoints)
        createMapFromBreakpoint(bp);
    mergeOverlappingMaps();
    for (VariantMap & m : this->refMaps)
        m.calculateMapLength();
    this->variantMap.calculateMapLength();
}

void VariantMapManager::fillMap(VariantMap & map, std::vector<Junction> & junctions)
{
    map.refName = this->chromosome;
    fillPositionsLeft(map, junctions);
    fillPositionsBetweenAllJunctions(map, junctions);
    fillPositionsRight(map, junctions);
}

void VariantMapManager::createMapFromBreakpoint(Breakpoint bp)
{
    VariantMap map;
    fillMap(map, bp);
    this->refMaps.push_back(map);
}

void VariantMapManager::fillMap(VariantMap & map, Breakpoint & bp)
{
    int startPosition = std::max(0, bp.getPosition() - 10*this->filterMargin);
    int endPosition = bp.getPosition() + 10*this->filterMargin;
    if (this->contigLengths.find(bp.getReferenceName()) != this->contigLengths.end()) 
        endPosition = std::min(endPosition, this->contigLengths[bp.getReferenceName()]);
    

    map.refName = bp.getReferenceName();
    map.regions.push_back(
        GenomicRegion(
            bp.getReferenceName(),
            startPosition,
            endPosition
        )
    );
}

void VariantMapManager::fillPositionsLeft(VariantMap & map, std::vector<Junction> & junctions)
{
    Junction startJunction = junctions[0];

    if (!startJunction.hasValidLeftSide())
        return;

    int startPosition = startJunction.getPositionLeft() + 10 * this->filterMargin * startJunction.getDirectionLeft();
    if (this->contigLengths.find(startJunction.getRefNameLeft()) != this->contigLengths.end())
        startPosition = std::max(std::min(startPosition, this->contigLengths[startJunction.getRefNameLeft()]), 0);
    else
        startPosition = std::max(startPosition, 0);
    bool reverse = startJunction.getDirectionLeft() > 0;

    GenomicRegion leftRegion(
        startJunction.getRefNameLeft(), 
        std::min(startPosition, startJunction.getPositionLeft()),
        std::max(startPosition, startJunction.getPositionLeft()),
        reverse
    );
    map.regions.push_back(leftRegion);
}

void VariantMapManager::fillPositionsBetweenAllJunctions(VariantMap & map, std::vector<Junction> & junctions)
{
    for (unsigned i = 0; i < junctions.size() - 1; ++i)
        fillPositionsBetweenJunctions(
            map, 
            junctions[i], 
            junctions[i+1], 
            i+1
        );
}

void VariantMapManager::fillPositionsBetweenJunctions(VariantMap & map, Junction & leftJunction, Junction & rightJunction, int sectionIndex)
{
    if (leftJunction.getDirectionRight() == rightJunction.getDirectionLeft() || leftJunction.getRefNameRight() != rightJunction.getRefNameLeft())
        throw std::runtime_error("Variant Disjointed. Cannot create mapping.");
    bool reverse = leftJunction.getDirectionRight() < 0;

    GenomicRegion region(
        leftJunction.getRefNameRight(),
        std::min(leftJunction.getPositionRight(), rightJunction.getPositionLeft()),
        std::max(leftJunction.getPositionRight(), rightJunction.getPositionLeft()),
        reverse
    );
    map.regions.push_back(region);
}

void VariantMapManager::fillPositionsRight(VariantMap & map, std::vector<Junction> & junctions)
{
    Junction endJunction = junctions[junctions.size()-1];
    if (!endJunction.hasValidRightSide())
        return;

    int endPosition = endJunction.getPositionRight() + 10*this->filterMargin*endJunction.getDirectionRight();
    if (this->contigLengths.find(endJunction.getRefNameRight()) != this->contigLengths.end())
        endPosition = std::max(0, std::min(endPosition, this->contigLengths[endJunction.getRefNameRight()]));
    else
        endPosition = std::max(0, endPosition);
    int reverse = endJunction.getDirectionRight() < 0;

    GenomicRegion rightRegion(
        endJunction.getRefNameRight(),
        std::min(endJunction.getPositionRight(), endPosition),
        std::max(endJunction.getPositionRight(), endPosition),
        reverse
    );
    map.regions.push_back(rightRegion);
}

void VariantMapManager::mergeOverlappingMaps()
{
    for (auto it = this->refMaps.begin(); it != this->refMaps.end();)
    {
        GenomicRegion r1 = it->regions[0];
        bool merged = false;
        for (auto it1 = it + 1; it1 != this->refMaps.end(); )
        {
            GenomicRegion r2 = it1->regions[0];
            if (r1.overlaps(r2))
            {
                it->regions[0].mergeWithRegion(r2);
                it1 = this->refMaps.erase(it1);
                merged = true;
                break;
            } else {
                ++it1;
            }
        }
        if (!merged)
            ++it;
    }

    if (this->variantMapPresent)
    {
        for (auto it = this->refMaps.begin(); it != this->refMaps.end(); )
        {
            GenomicRegion r1 = it->regions[0];
            if (r1.overlaps(this->refRegionLeft) || r1.overlaps(this->refRegionRight))
            {
                bool toErase = false;
                if (r1.overlaps(this->refRegionLeft))
                {
                    if (r1.overlaps(this->variantMap.regions[0]))
                    {
                        if (this->variantMap.regions[0].isReverse())
                            this->variantMap.regions[0].setRegionEnd(std::max(this->variantMap.regions[0].getRegionEnd(), r1.getRegionEnd()));
                        else
                            this->variantMap.regions[0].setRegionStart(std::min(this->variantMap.regions[0].getRegionStart(), r1.getRegionStart()));
                        createReferenceRegions();
                        toErase = true;
                    }
                }
                if (r1.overlaps(this->refRegionRight))
                {
                    int regionIdx = this->variantMap.regions.size() - 1;
                    if (r1.overlaps(this->variantMap.regions[regionIdx]))
                    {
                        if (! this->variantMap.regions[regionIdx].isReverse())
                            this->variantMap.regions[regionIdx].setRegionEnd(std::max(this->variantMap.regions[regionIdx].getRegionEnd(), r1.getRegionEnd()));
                        else
                            this->variantMap.regions[regionIdx].setRegionStart(std::min(this->variantMap.regions[regionIdx].getRegionStart(), r1.getRegionStart()));
                        createReferenceRegions();
                        toErase = true;
                    }
                }
                if (!toErase)
                    ++it;
                else
                    it = this->refMaps.erase(it);
            } else {
                it = this->refMaps.erase(it);
            }
        }
    }
}

void VariantMapManager::createReferenceRegions()
{
    std::string rLeft = this->variantMap.regions[0].getReferenceName();
    int dLeft = (this->variantMap.regions[0].isReverse() ? 1 : -1);
    int firstPosition = (this->variantMap.regions[0].isReverse() ? this->variantMap.regions[0].getRegionEnd() : this->variantMap.regions[0].getRegionStart());
    int leftRegionEnd = firstPosition + dLeft; //* (this->filterMargin * 10 + 1);

    int lastIdx = this->variantMap.regions.size() - 1;
    std::string rRight = this->variantMap.regions[lastIdx].getReferenceName();
    int dRight = this->variantMap.regions[lastIdx].isReverse() ? -1 : 1;
    int lastPosition = this->variantMap.regions[lastIdx].isReverse() ? this->variantMap.regions[lastIdx].getRegionStart() : this->variantMap.regions[lastIdx].getRegionEnd();
    int rightRegionEnd = lastPosition + dRight;// * (this->filterMargin * 10 + 1);

    if (dLeft < 0) {
        this->refRegionLeft = GenomicRegion(rLeft, 0, leftRegionEnd);
    } else {
        if (this->contigLengths.find(rLeft) == this->contigLengths.end())
        {
            this->refRegionLeft = GenomicRegion(rLeft, leftRegionEnd, 300000000);
        } else {
            this->refRegionLeft = GenomicRegion(rLeft, leftRegionEnd, this->contigLengths[rLeft]);
        }
    }

    if (dRight > 0) {
        if (this->contigLengths.find(rRight) == this->contigLengths.end())
            this->refRegionRight = GenomicRegion(rRight, rightRegionEnd, 300000000);
        else
            this->refRegionRight = GenomicRegion(rRight, rightRegionEnd, this->contigLengths[rRight]);
    } else {
        this->refRegionRight = GenomicRegion(rRight, 0, rightRegionEnd);
    }
}

std::vector<VariantMap> VariantMapManager::getMaps()
{
    std::vector<VariantMap> maps;
    if (this->variantMapPresent)
        maps.push_back(this->variantMap);
    for (auto m : this->refMaps)
        maps.push_back(m);
    return maps;
}

void VariantMapManager::print()
{
    for (auto m : getMaps())
    {
        for (GenomicRegion region : m.regions)
            region.print();
        std::cout << std::endl;
    }
}

seqan::String<seqan::Dna5String> VariantMapManager::getChromosomeSequences(SeqFileHandler & seqFileHandler)
{
    seqan::String<seqan::Dna5String> sequences;
    if (this->variantMapPresent)
    {
        seqan::Dna5String variantSequence;
        for (GenomicRegion & region : this->variantMap.regions)
        {
            region.readSequence(seqFileHandler);
            seqan::append(variantSequence, region.getSequence());
        }
        seqan::appendValue(sequences, variantSequence);
    }
    for (VariantMap & map : this->refMaps)
    {
        seqan::Dna5String seq;
        for (GenomicRegion & region : map.regions)
        {
            region.readSequence(seqFileHandler);
            seqan::append(seq, region.getSequence());
        }
        seqan::appendValue(sequences, seq);
    }
    return sequences;
}