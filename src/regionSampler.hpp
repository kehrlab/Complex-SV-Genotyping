#ifndef REGIONSAMPLERHEADER
#define REGIONSAMPLERHEADER

#include <vector>
#include <random>

#include "custom_types.hpp"
#include "genomicRegion.hpp"
#include "bamFileHandler.hpp"
#include "genomicRegion.hpp"

class RegionSampler
{
    std::vector<ContigInfo> contigInfo;
    std::vector<GenomicRegion> insertSizeRegions;
    std::vector<GenomicRegion> gcBiasRegions;

    std::vector<std::string> commonContigs;
    std::vector<int> minContigLengths;

    uint32_t intervalSizes;
    uint32_t totalIntervalLength;
    bool gcBias;

    void determineCommonContigs();
    std::vector<GenomicRegion> createGenomicRegions();
    void subsampleRegions(std::vector<GenomicRegion>);
    bool allRegionsValid(std::vector<GenomicRegion>);
    void sampleInsertSizeRegions(std::vector<GenomicRegion>);

    public:
    RegionSampler();
    RegionSampler(std::vector<ContigInfo>, std::vector<GenomicRegion>);
    std::vector<GenomicRegion> & sampledInsertRegions();
};

#endif