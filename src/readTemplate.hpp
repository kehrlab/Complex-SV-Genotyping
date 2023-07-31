#ifndef READTEMPLATEHEADER
#define READTEMPLATEHEADER

#include <seqan/bam_io.h>
#include "bamFileHandler.hpp"
#include "record.hpp"
#include "filter.hpp"
#include "junction.hpp"
#include <vector>
#include <unordered_set>

class ReadTemplate
{
    std::vector<BamRecord> records;
    std::vector<int> firstIndices;
    std::vector<int> lastIndices;
    int primaryFirst, primaryLast;
    int fivePrimeFirst, fivePrimeLast;
    std::string templateName;

    std::string orientation;
    int insertSize;
    bool split;
    bool splitFirst;
    bool splitLast;
    bool spanning;
    bool interChromosome;

    bool suspectedSplit;

    float templateGCContent;

    float templateWeight;
    std::vector<std::vector<std::string>> interChromosomeNames;
    std::unordered_set<int> overlappingRegions;
    std::string regionString;
    std::string junctionString;
    std::string bpString;
    std::unordered_set<int> spanningBreakpoints;
    std::unordered_set<int> splittingJunctions;

    void determineInformativeRecords();
    void calculateInsertSize();
    void determineOrientation();
    
    int softClippingLeftLength(BamRecord &);
    int softClippingRightLength(BamRecord &);
    bool isClippedAtJunction(int &, int &, BamRecord &, Junction &);
    void findSplitsBasedOnClipping(std::vector<Junction> &);
    void findSplitsBasedOnGaps(std::vector<Junction> &);
    void findClippedSplitsOnChromosome(std::string, std::vector<Junction> &);
    bool alignsWithinExpectedDistance(std::vector<Junction> &, int, BamRecord &, BamRecord &);
    bool alignsWithinRegion(int &, GenomicRegion &, BamRecord &, int);
    bool readSpansJunction(BamRecord &, Junction &);
    
    public:
    ReadTemplate();
    ReadTemplate(std::vector<BamRecord>);
    ReadTemplate(std::vector<BamRecord>, std::vector<Junction> &, std::vector<Breakpoint> &);

    void determineGroup(std::vector<Junction> &, std::vector<Breakpoint> &, std::vector<GenomicRegion> &);
    void findSplitReads(std::vector<Junction> &);
    void findSpanningReads(std::vector<Breakpoint> &);
    void determineOverlappingRegions(std::vector<GenomicRegion> &);
    void determineLocationStrings();
    void determineFivePrimeEnds();
    void markSuspectedSplit();



    std::string getOrientation();
    int getInsertSize();
    bool containsSplitRead();
    bool containsSpanningRead();
    bool alignsAcrossChromosomes();

    std::string getPrimaryRecordString();
    bool isProperPair();
    bool isInterestingReadPair(ReadPairFilter &);

    std::vector<BamRecord> & getRecords();
    std::vector<BamRecord> getPrimaryRecords();
    float getTemplateWeight();
    std::vector<std::vector<std::string>> getInterChromosomeNames();
    void print();
    void setGCContent(float);
    float getGCContent();
    bool containsSuspectedSplit();
    std::string getName();
    std::string getRegionString();
    std::string getJunctionString();
    std::string getBreakpointString();
    std::vector<int> getMappingQualities();
};

#endif
