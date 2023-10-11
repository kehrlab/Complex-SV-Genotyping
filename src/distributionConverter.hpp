#ifndef DISTRIBUTIONCONVERTERHEADER
#define DISTRIBUTIONCONVERTERHEADER

#include <vector>
#include <string>
#include "custom_types.hpp"
#include "filter.hpp"
#include "genotypeDistribution.hpp"
#include "libraryDistribution.hpp"
#include "options.hpp"
#include "variant.hpp"
#include "readTemplate.hpp"
#include <random>
#include "variantProfile.hpp"

class DistributionConverter
{
    ReadPairFilter filter;
    LibraryDistribution sampleDistribution;
    complexVariant variant;
    std::string filename;

    ProgramOptions options;

    int readLength;
    float lowerInsertLimit;
    float upperInsertLimit;

    std::vector<GenotypeDistribution> genotypeDistributions;
    std::vector<GenotypeDistribution> alleleDistributions;
    std::vector<std::string> genotypeNames;
    std::vector<std::string> alleleNames;
    std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> chromosomeStructures;

    void addVariantRegions(std::string, GenotypeDistribution &, Allele &, bool);
    void addVariantRegionsWithGCCorrection(std::string, GenotypeDistribution &, Allele &);
    void addVariantRegionsWithoutGCCorrection(std::string, GenotypeDistribution &, Allele &);
    float determineMinProbability();
    void setMinProbabilities(float);

    float extractGCCorrectionFactor(float, int);
    void determineDifficulty();

    public:
    DistributionConverter();
    DistributionConverter(complexVariant, std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> &, LibraryDistribution, int, BamFileHandler &, ProgramOptions &);
    DistributionConverter(VariantProfile &, LibraryDistribution &, BamFileHandler &, bool);

    void initDistributions(BamFileHandler &, int);
    void createDistributions();
    void createVariantDistributions();
    void createAlleleDistribution(GenotypeDistribution &, Allele &);
    void createMixedDistributions();
    void scaleDistributions();
    int addSimulatedTemplateToDistribution(VariantMap &, int, int, float, GenotypeDistribution &, Allele &);

    ReadPairFilter & getReadPairFilter();

    std::vector<GenotypeDistribution> & getGenotypeDistributions();
    std::vector<std::string> & getGenotypeNames();
};

#endif
