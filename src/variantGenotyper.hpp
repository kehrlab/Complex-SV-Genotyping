#ifndef VARIANTGENOTYPEHEADER
#define VARIANTGENOTYPEHEADER

#include "genotypeDistribution.hpp"
#include "genotypeResult.hpp"
#include "likelihoodCalculator.hpp"
#include "recordManager.hpp"
#include "variant.hpp"
#include "variantProfile.hpp"
#include <boost/filesystem.hpp>

class VariantGenotyper
{
    // references
    LibraryDistribution & sampleDistribution;
    complexVariant & variant;
    ProgramOptions & options;

    // sample-specific resources
    std::string bamFileName;
    BamFileHandler bamFile;
    RecordManager bamRecords;
    
    // result
    GenotypeResult variantGenotype;

    void createGenotypeDistribution();
    void calculateGenotypeLikelihoods();
    void writeInsertSizeDistributions(LikelihoodCalculator &);

    public:
    VariantGenotyper(complexVariant &, std::string, LibraryDistribution &, ProgramOptions &);
    void genotype();
    void genotype(VariantProfile &);
    GenotypeResult getResult();
    std::string getVariantName();
    std::string getSampleName();
};

#endif
