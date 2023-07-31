#ifndef PROGRAMOPTIONSHEADER
#define PROGRAMOPTIONSHEADER

#include "genomicRegion.hpp"
#include <string>
#include <vector>

class ProgramOptions
{
    std::string outFile;
    std::string variantFile;
    std::string vcfFile;
    std::string refFile;

    std::string sequenceDirectory;

    std::vector<std::string> bamFileNames;
    std::vector<GenomicRegion> samplingRegions;

    bool verbose;
    bool output;
    bool wholeGenome;
    bool profile;
    int nThreads;
    int minQ;
    bool outputDistributions;
    int estimateDiffQ;
    int mode;
    bool splitReads;
    bool spanningReads;
    bool normalReads;
    bool noInsertSizes;
    bool useQualities;
    bool gcCorrect;
    bool loadToMemory;
    int Ts, Tv;
    int simCoverage;
    bool stats;

    void determineBamFiles(std::string, std::string);
    void determineSamplingRegions(std::string);

    public:
    ProgramOptions();
    ProgramOptions(
        std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string,
        bool, int, int, int, int, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, int, bool
    );
    std::vector<std::string> getBamFileNames();
    std::vector<GenomicRegion> getSamplingRegions();
    std::string getOutputFileName();
    std::string getVariantFileName();
    std::string getVCFName();
    std::string getRefFileName();
    std::string getSequenceDirectory();

    bool isOptionVerbose();
    bool isOptionWriteUsedReads();
    bool isOptionUseWholeGenome();
    bool isOptionWriteToFile();
    bool isOptionOutputDistributions();
    int getEstimateDiffQual();
    int getNumberOfThreads();
    int getMinimumMappingQuality();
    bool isOptionNoSplitReads();
    bool isOptionNoSpanningReads();
    bool isOptionNoNormalReads();
    bool isOptionProfile();
    bool isOptionNoInsertSizes();
    bool isOptionUseQualities();
    int getDistributionMode();
    bool isOptionGCCorrect();
    bool isOptionLoadToMemory();
    void determineNumThreads(int, int);
    int getNumberOfSampleThreads();
    int getNumberOfVariantThreads();
    bool isOptionSample();
    int getCoverage();
    bool isOptionStats();
};

#endif
