#ifndef VARIANTMAPHEADER
#define VARIANTMAPHEADER

#include <string>
#include "genomicRegion.hpp"
#include "readTemplate.hpp"
#include "record.hpp"

class VariantMap
{
    public:
    VariantMap();
    
    std::string refName;
    std::vector<GenomicRegion> regions;
    std::vector<int> regionLengths;
    std::vector<int> regionIndices;
    
    unsigned totalLength;
    bool possibleSplitRead;

    void calculateMapLength();
    ReadTemplate simulateTemplate(int, int, int);
    void getGCContentVector(std::vector<float> &, unsigned, int, int);
    seqan::Dna5String getBaseAtIndex(int);
    GenomicRegion getPositionOnReference(int);
    int getRegionIndex(int);
    std::vector<BamRecord> getReadsFromIndices(int, int, bool, bool);
};

#endif
