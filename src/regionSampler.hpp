#ifndef REGIONSAMPLERHEADER
#define REGIONSAMPLERHEADER

#include "custom_types.hpp"
#include <vector>
#include "genomicRegion.hpp"
#include "options.hpp"

class RegionSampler
{
    std::vector<ContigInfo> contigInfo;
    std::vector<GenomicRegion> insertSizeRegions;
    std::vector<GenomicRegion> gcBiasRegions;

    std::vector<std::string> commonContigs;
    std::vector<int> minContigLengths;

    int intervalSizes;
    int totalIntervalLength;
    bool gcBias;

    void determineCommonContigs();
    std::vector<GenomicRegion> createGenomicRegions();
    void subsampleRegions(std::vector<GenomicRegion>);
    bool allRegionsValid(std::vector<GenomicRegion>);
    void sampleInsertSizeRegions(std::vector<GenomicRegion>);

    public:
    RegionSampler();
    RegionSampler(std::vector<ContigInfo>, std::vector<GenomicRegion>);
    RegionSampler(std::vector<ContigInfo>, ProgramOptions &, std::vector<GenomicRegion>);
    std::vector<GenomicRegion> & sampledInsertRegions();
};

#endif