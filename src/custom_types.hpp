#ifndef CUSTOM_TYPE_HEADER
#define CUSTOM_TYPE_HEADER

#include <seqan/sequence.h> 
#include <seqan/bam_io.h>
#include <string>
#include <vector>

#include "genomicRegion.hpp"
#include "junction.hpp"
#include "breakpoint.hpp"

struct ContigInfo
{
    std::vector<std::string> cNames;
    std::vector<int> cLengths;
};

struct TemplatePosition
{
    int begin, end;
    std::string chr;
};

struct JunctionRegion
{
    std::vector<GenomicRegion> regions;
    std::vector<int> junctionIndices;
    std::vector<int> breakpointIndices;
    std::vector<Junction> junctions;
    std::vector<Breakpoint> breakpoints;
    int length;
    std::string chromosome;
};

struct BreakpointRegion
{
    GenomicRegion region;
    std::vector<Breakpoint> breakpoints;
};

struct VariantRegions
{
    std::vector<JunctionRegion> regions;
    std::vector<int> distanceFromLast;
    std::vector<int> distanceToNext;
};

struct SplitAlignmentInfo
{
    std::vector<std::vector<int>> junctionIndices;
    std::vector<int> insertSize; 
};

#endif
