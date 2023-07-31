#ifndef PARSERHEADER
#define PARSERHEADER

#define SEQAN_ENABLE_PARALLELISM 1
#include "custom_types.hpp"
#include "seqan/arg_parse.h"
#include "options.hpp"

class parser {
    seqan::ArgumentParser argParser;
    seqan::ArgumentParser::ParseResult parsingStatus;
    ProgramOptions options;

    void addOptionsIO();
    void addOptionsParams();

    void addOptionInputFile();
    void addOptionInputFileList();
    void addOptionOutputFile();
    void addOptionConfigFile();
    void addOptionOutputVCF();
    void addOptionGenomeFile();
    void addOptionSequenceDirectory();
    void addOptionSamplingRegions();

    void addOptionVerbose();
    void addOptionNumThreads();
    void addOptionWholeGenome();
    void addOptionProfiling();
    void addOptionEstimateDifficulty();
    void addOptionNoSplits();
    void addOptionNoSpanning();
    void addOptionNoStandard();
    void addOptionNoInsertSize();
    void addOptionUseQualities();
    void addOptionMinQuality();
    void addOptionMode();
    void addOptionGCCorrect();
    void addOptionLoadToMemory();
    void addOptionCoverage();
    void addOptionStats();

    void addOptionOutputInsertSizeDistributions();

    void parseOptions(int, char const **);
    
    seqan::ArgumentParser::ParseResult getArguments();

    void checkForArgumentConflict();
    bool singleFileAndListPresent();
    bool singleFileAndListMissing();
    void checkParsingStatus();
    void extractOptions();
    void checkGCRequirements();

    public:
    parser();
    parser(int, char const **);
    bool wasSuccessful();
    ProgramOptions getOptions();
};

#endif
