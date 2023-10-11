#ifndef GENOTYPERHEADER
#define GENOTYPERHEADER

#include "genotypeResult.hpp"
#include "libraryDistribution.hpp"
#include "parser.hpp"
#include "regionSampler.hpp"
#include "variant.hpp"
#include "options.hpp"
#include "bamFileHandler.hpp"
#include "genotypeResult.hpp"
#include <unordered_map>
#include <fstream>
#include <chrono>
#include <boost/filesystem.hpp>
#include "sample.hpp"
#include "variantGenotyper.hpp"
#include "vcfWriter.hpp"


class Genotyper
{
    ProgramOptions options;
    std::vector<complexVariant> variants;
    RegionSampler regionSampler;
    SeqFileHandler referenceFile;

    std::vector<std::string> fileNames;
    std::vector<std::vector<GenotypeResult>> genotypeResults;
    
    std::vector<LibraryDistribution> sampleDistributions;
    std::vector<std::string> sampleNames;

    std::vector<VariantProfile> variantProfiles;

    std::vector<ContigInfo> contigInfos;

    int maxFilterMargin;
    int maxReadLength;
    std::unordered_map<std::string, int> contigLengths;

    std::ofstream outputFile;

    bool writeToFile = false;

    void createSamples();
    void createResultVector();

    public:
    Genotyper();
    Genotyper(ProgramOptions);
    void setOptions(ProgramOptions &);
    void getBamFileNames();
    void createVariants();
    void printVariantDetails();
    void extractContigInfos();
    void sampleRegionsForInsertSizeDistribution();
    void genotypeAllSamples();
    void writeResults();
    void writeToVCF();
    void writeInsertSizeDistributions(complexVariant &, LikelihoodCalculator &);
    void openReferenceFile();
    void writeVariantAlleles();
    void closeReferenceFile();
    void createAlleleMaps();
    void gatherContigLengths();
    void writeStats();
    void createVariantProfiles();
};

#endif
