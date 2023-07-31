#include "filter.hpp"

ReadPairFilter::ReadPairFilter(std::vector<Breakpoint> const & breakpoints, float mean, float sd)
{
    this->margin = calculateBreakpointMargin(mean, sd);
    defineBreakpointRegions(breakpoints);
}

int ReadPairFilter::calculateBreakpointMargin(float mean, float sd)
{
    return (int)(mean + 2*sd);
}

void ReadPairFilter::defineBreakpointRegions(const std::vector<Breakpoint> & breakpoints)
{
    for (auto bp : breakpoints)
        addBreakpointRegion(bp);
    mergeBreakpointRegions();
}

void ReadPairFilter::addBreakpointRegion(Breakpoint breakpoint)
{
    GenomicRegion bpRegion(
        breakpoint.getReferenceName(), 
        breakpoint.getPosition() - this->margin, 
        breakpoint.getPosition() + this->margin
    );
    this->variantRegions.push_back(bpRegion);
}

void ReadPairFilter::mergeBreakpointRegions()
{
    for (auto it = this->variantRegions.begin(); it != this->variantRegions.end();)
    {
        bool merged = false;
        for (auto it1 = it + 1; it1 != this->variantRegions.end(); )
        {
            if (it1->overlaps(*it))
            {
                it->mergeWithRegion(*it1);
                it1 = this->variantRegions.erase(it1);
                merged = true;
                break;
            } else {
                ++it1;
            }
        }
        if (!merged)
            ++it;
    }
}

bool ReadPairFilter::fragmentPassesFilter(std::string rName, int begin, int end)
{
    GenomicRegion pairRegion(rName, begin, end);
    for (auto region : this->variantRegions)
        if (region.overlaps(pairRegion))
            return true;
    return false;
}

bool ReadPairFilter::fragmentPassesFilter(BamRecord & firstRecord, BamRecord & lastRecord)
{
    GenomicRegion firstRegion = firstRecord.getAlignmentRegion();
    GenomicRegion lastRegion = lastRecord.getAlignmentRegion();
    for (auto region : this->variantRegions)
        if (region.overlaps(firstRegion) || region.overlaps(lastRegion))
            return true;

    if (firstRecord.getReferenceName() == lastRecord.getReferenceName())
    {
        if (fragmentPassesFilter(firstRecord.getReferenceName(), std::min(firstRegion.getRegionStart(), lastRegion.getRegionStart()), std::max(firstRegion.getRegionEnd(), lastRegion.getRegionEnd())))
            return true;
    }

    return false;
}

int ReadPairFilter::getMargin()
{
    return this->margin;
}