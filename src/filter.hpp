#ifndef READPAIRFILTERHEADER
#define READPAIRFILTERHEADER

#include <vector>
#include "breakpoint.hpp"
#include "genomicRegion.hpp"
#include "record.hpp"

class ReadPairFilter
{
    std::vector<GenomicRegion> variantRegions;
    int margin;

    public:
    ReadPairFilter() : variantRegions(){};
    ReadPairFilter(std::vector<Breakpoint> const &, float, float);
    static int calculateBreakpointMargin(float, float);
    void defineBreakpointRegions(std::vector<Breakpoint> const &);
    void addBreakpointRegion(Breakpoint);
    void mergeBreakpointRegions();
    bool fragmentPassesFilter(std::string, int, int);
    bool fragmentPassesFilter(BamRecord &, BamRecord &);
    int getMargin();
    bool fragmentOverlapsNovelJunction(BamRecord &, BamRecord &);
};

#endif