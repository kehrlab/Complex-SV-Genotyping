#ifndef READTEMPLATEHEADER
#define READTEMPLATEHEADER

#include <seqan/bam_io.h>
#include <vector>
#include <unordered_set>
#include <string>

#include "bamFileHandler.hpp"
#include "record.hpp"
#include "filter.hpp"
#include "junction.hpp"
#include "custom_types.hpp"
#include "genomicRegion.hpp"

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
    std::unordered_set<int> overlappingRegions;
    std::string junctionString;
    std::string chromosomeString;
    std::string bpString;
    std::unordered_set<int> spanningBreakpoints;
    std::unordered_set<int> splittingJunctions;

    void determineInformativeRecords();
    void calculateInsertSize();
    void determineOrientation();
    
    int softClippingLeftLength(BamRecord &);
    int softClippingRightLength(BamRecord &);
    void findSplitsBasedOnClipping(std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> &);
    void findSplitsBasedOnGaps(std::vector<Junction> &);
    void findClippedSplitsOnChromosome(std::string, JunctionRegion &, SplitAlignmentInfo &);
    bool alignsWithinExpectedDistance(std::vector<Junction> &, uint32_t, BamRecord &, BamRecord &);
    bool alignsWithinRegion(int &, GenomicRegion &, BamRecord &, int);
    bool readSpansJunction(BamRecord &, Junction &);

    void getIndexRegions(std::vector<GenomicRegion> &, std::vector<std::vector<int>> &, BamRecord &, JunctionRegion &);
    
    public:
    ReadTemplate();
    ReadTemplate(std::vector<BamRecord>);

    void findSplitReads(std::vector<Junction> &, std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> &);
    void findSpanningReads(std::vector<Breakpoint> &);
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
    void print();
    void setGCContent(float);
    float getGCContent();
    bool containsSuspectedSplit();
    std::string getName();
    std::string getJunctionString();
    std::string getBreakpointString();
    std::string getChromosomeString();
    std::vector<int> getMappingQualities();
};

#endif
