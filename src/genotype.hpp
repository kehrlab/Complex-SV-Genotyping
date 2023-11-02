#ifndef GENOTYPEHEAEDER
#define GENOTYPEHEAEDER

#include <seqan/arg_parse.h>
#include <stdexcept>

#include "filter.hpp"
#include "genotypeDistribution.hpp"
#include "genotypeResult.hpp"
#include "libraryDistribution.hpp"
#include "recordManager.hpp"
#include "sample.hpp"
#include "variantProfile.hpp"


struct genotypeParameters
{
    std::string variantList;
    std::string sampleList;

    std::string outputPrefix;
    std::string variantFile;
    std::string vcfFile;

    std::string priorFile;

    int nThreads, minMapQ;
    bool difficulties, distributions;

    genotypeParameters():
    variantList(""), sampleList(""), outputPrefix(""), variantFile(""),
    vcfFile(""), priorFile(""), nThreads(1), minMapQ(0), difficulties(false), 
    distributions(false)
    {}
};

int genotype(int argc, const char ** argv);

seqan::ArgumentParser::ParseResult parseGenotypeArgs(seqan::ArgumentParser & argParser, int argc, const char ** argv);

genotypeParameters getGenotypeParameters(seqan::ArgumentParser & argParser);

inline void loadVariantProfiles(std::vector<VariantProfile> &, genotypeParameters &, std::vector<std::string> &, int);
inline void checkProfileParameters(int &, int &, int &, std::vector<VariantProfile> &, genotypeParameters &);
inline void getGenotypePriors(std::unordered_map<std::string, std::unordered_map<std::string, float>> & genotypePriors, genotypeParameters & params);
inline void checkSampleParameters(std::vector<std::string> &, int &, int &, int &, genotypeParameters &);
inline void loadReadPairs(VariantProfile &, RecordManager &, BamFileHandler &, Sample &);
inline void createGenotypeDistributions(
    VariantProfile &, std::vector<std::string> &, std::vector<GenotypeDistribution> &, 
    float eps, std::vector<std::string> &, LibraryDistribution &);
inline void calculateGenotypeLikelihoods(
    VariantProfile & variantProfile, std::vector<std::string> & genotypeNames, 
    std::vector<GenotypeDistribution> & genotypeDistributions, 
    RecordManager & bamRecords, ReadPairFilter & filter, GenotypeResult & result
    );

inline void adjustLikelihoods(ReadTemplate & rT, std::vector<std::string> & names, std::vector<GenotypeDistribution> & distributions, GenotypeResult & result);

void writeInsertSizeDistributions(complexVariant & variant, GenotypeResult & result, std::vector<std::string> & names, std::vector<GenotypeDistribution> & distributions);

void calculateDifficulty(float &, std::vector<std::string> & genotypeNames, std::vector<GenotypeDistribution> & genotypeDistributions);


#endif