#ifndef GENOMICREGIONHEADER
#define GENOMICREGIONHEADER

#include <string>
#include <iostream>
#include <vector>

#include "seqFileHandler.hpp"


class GenomicRegion
{
    std::string rName;
    int64_t begin, end;

    bool reverse;
    bool containsSequence;
    seqan::Dna5String sequence;

    float nFraction;
    int longestN;
    int nStart, nEnd;

    void setReverse(bool);
    void determineNStats(SeqFileHandler &);

    public:
    GenomicRegion();
    GenomicRegion(std::string);
    GenomicRegion(std::string, int64_t, int64_t);
    GenomicRegion(std::string, int64_t, int64_t, bool);
    bool operator==(GenomicRegion);
    std::string getReferenceName();
    int64_t getRegionStart();
    int64_t getRegionEnd();
    void setReferenceName(std::string);
    void setRegionStart(int64_t);
    void setRegionEnd(int64_t);
    bool overlaps(GenomicRegion);
    bool overlaps(GenomicRegion, int);
    void mergeWithRegion(GenomicRegion);
    bool sequenceInMemory();
    void reverseRegion();
    void print();
    void clearSequence();
    bool isReverse();
    void readSequence(SeqFileHandler &);
    static std::vector<GenomicRegion> joinRegions(std::vector<GenomicRegion>);
    seqan::Dna5String & getSequence();
    bool isValidSample(SeqFileHandler &);
    std::string getRegionString();
};

#endif
