#ifndef LIKELIHOODCALCULATORHEADER
#define LIKELIHOODCALCULATORHEADER

#include "bamFileHandler.hpp"
#include "breakpoint.hpp"
#include "distributionConverter.hpp"
#include "genotypeDistribution.hpp"
#include "genotypeResult.hpp"
#include "junction.hpp"
#include "libraryDistribution.hpp"
#include "options.hpp"
#include "readTemplate.hpp"
#include "recordManager.hpp"
#include "variant.hpp"
#include "variantProfile.hpp"
#include <unordered_map>
#include <vector>

class LikelihoodCalculator
{
    std::vector<ReadTemplate> templates;
    complexVariant & variant;
    LibraryDistribution & sampleDistribution;
    BamFileHandler & bamFileHandler;
    ProgramOptions & options;
    int maxReadLength;

    float splitReadFactor;
    float spanningReadFactor;
    float normalReadFactor;
    bool useInsertSizes;

    float lowerInsertLimit;
    float upperInsertLimit;

    GenotypeResult result;

    std::vector<GenotypeDistribution> genotypeDistributions;
    std::vector<std::string> genotypeNames;
    ReadPairFilter filter;

    public:
    LikelihoodCalculator(RecordManager &, BamFileHandler &, complexVariant &, LibraryDistribution &, ProgramOptions &);
    void createInsertSizeDistributions(std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> &);
    void createInsertSizeDistributions(VariantProfile &);
    void calculateLikelihoods(std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> &);
    void calculateLikelihoods(VariantProfile &);
    void adjustLikelihoods(ReadTemplate &);
    std::vector<GenotypeDistribution> & getDistributions();
    std::vector<std::string> & getGenotypeNames();
    void clearData();

    GenotypeResult getResult();
};

#endif
