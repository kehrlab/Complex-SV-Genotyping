#ifndef PROFILESAMPLESHEADER
#define PROFILESAMPLESHEADER

#include "genomicRegion.hpp"
#include <seqan/arg_parse.h>
#include "bamFileHandler.hpp"
#include "custom_types.hpp"
#include "regionSampler.hpp"
#include "sample.hpp"
#include "seqFileHandler.hpp"
#include <fstream>
#include <stdexcept>

struct sampleProfileParams {
    std::vector<std::string> sampleFiles;
    std::string outFile;
    std::string outDir;
    std::vector<GenomicRegion> regions;

    bool wholeGenome;
    int nThreads;

    sampleProfileParams ():
    sampleFiles(std::vector<std::string>()),
    outFile(""),
    outDir(""),
    regions(std::vector<GenomicRegion>()),
    wholeGenome(false),
    nThreads(1)
    {}
};

int profileSamples(int argc, const char ** argv);

seqan::ArgumentParser::ParseResult parseSampleProfileArgs(seqan::ArgumentParser & argParser, int argc, const char ** argv);
sampleProfileParams getSampleProfileParameters(const seqan::ArgumentParser & argParser);

#endif