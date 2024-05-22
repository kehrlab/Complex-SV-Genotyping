#ifndef PROFILEVARIANTSHEADER
#define PROFILEVARIANTSHEADER

#include <stdexcept>
#include <seqan/arg_parse.h>

#include "profileSamples.hpp"
#include "variantProfile.hpp"
#include "variantParser.hpp"
#include "variant.hpp"
#include "profileHandler.hpp"

struct variantProfileParams {
    std::vector<std::string> variantFileNames;
    std::string outFile;
    std::string outDir;

    std::string sampleProfileFile;
    int margin;
    int nThreads;

    bool popdel;

    variantProfileParams():
    variantFileNames (std::vector<std::string>()), outFile(""), outDir(""), sampleProfileFile (""),
    margin(500), nThreads(1), popdel(false)
    {} 
};

int profileVariants(int argc, const char ** argv);
inline void createVariantProfile(const complexVariant & variant, const variantProfileParams & params, int overlap, int readLength, int sMin, int sMax, ContigInfo cInfo, std::atomic_int64_t & counter, std::vector<bool> & availability, int idx);
inline void trackProgress(const std::atomic_int64_t & counter, int64_t targetSize);
seqan::ArgumentParser::ParseResult parseVariantProfileArgs(seqan::ArgumentParser & argParser, int argc, const char ** argv);
variantProfileParams getVariantProfileParameters(const seqan::ArgumentParser & argParser);
ContigInfo mergeContigLengths(const std::vector<std::unordered_map<std::string, int>> &);

#endif