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
    std::unordered_map<int, int> deletionIndices;
    
    unsigned totalLength;
    bool possibleSplitRead;

    void calculateMapLength();
    void detectSmallDeletions();
    
    ReadTemplate simulateTemplate(int, int, int);
    void getGCContentVector(std::vector<float> &, unsigned, int, int);

    private:
    bool isDeletion(int &, GenomicRegion &, GenomicRegion &);
    int getRegionIndex(int);
    std::vector<BamRecord> getReadsFromIndices(int, int, bool, bool);
    GenomicRegion getPositionOnReference(int);
    seqan::Dna5String getBaseAtIndex(int);
};

#endif
