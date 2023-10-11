#ifndef GENOMICREGIONHEADER
#define GENOMICREGIONHEADER

#include <string>
#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>

#include "seqFileHandler.hpp"


class GenomicRegion
{
    std::string rName;
    int begin, end;

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
    GenomicRegion(std::string, int, int);
    GenomicRegion(std::string, int, int, bool);
    bool operator==(GenomicRegion);
    std::string getReferenceName();
    int getRegionStart();
    int getRegionEnd();
    void setReferenceName(std::string);
    void setRegionStart(int);
    void setRegionEnd(int);
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
