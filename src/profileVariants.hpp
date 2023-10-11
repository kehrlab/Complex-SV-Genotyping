#ifndef PROFILEVARIANTSHEADER
#define PROFILEVARIANTSHEADER

#include <stdexcept>
#include <seqan/arg_parse.h>

#include "profileSamples.hpp"
#include "variantProfile.hpp"
#include "variantParser.hpp"
#include "variant.hpp"

struct variantProfileParams {
    std::string variantFile;
    std::string outFile;
    std::string outDir;

    std::string sampleProfileFile;
    int sMin;
    int sMax;
    int readLength;
    int margin;
    int nThreads;

    variantProfileParams():
    variantFile (""), outFile(""), outDir(""), sampleProfileFile (""),
    sMin(1), sMax(1000), readLength(100), margin(500)
    {} 
};

int profileVariants(int argc, const char ** argv);
seqan::ArgumentParser::ParseResult parseVariantProfileArgs(seqan::ArgumentParser & argParser, int argc, const char ** argv);
variantProfileParams getVariantProfileParameters(const seqan::ArgumentParser & argParser);
std::unordered_map<std::string, int> mergeContigLengths(const std::vector<std::unordered_map<std::string, int>> &);

#endif