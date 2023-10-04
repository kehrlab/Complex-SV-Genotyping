#ifndef SAMPLEHEADER
#define SAMPLEHEADER

#include <vector>
#include "bamFileHandler.hpp"
#include "genotypeDistribution.hpp"
#include "genotypeResult.hpp"
#include "recordManager.hpp"
#include "regionSampler.hpp"
#include "seqFileHandler.hpp"
#include "variant.hpp"
#include "variantGenotyper.hpp"
#include "libraryDistribution.hpp"
#include <chrono>

class Sample
{
    // global resources
    SeqFileHandler & referenceFile;
    ProgramOptions & options;
    std::vector<complexVariant> & variants;
    RegionSampler & regionSampler;

    // local resources
    std::string filename;
    std::string distributionDirectory;
    BamFileHandler bamFile;

    // results
    LibraryDistribution sampleDistribution;
    
    bool bamFileOpen;
    int maxReadLength;

    // private functions
    void calculateDefaultDistributions();
    void createSampleDistribution(std::unordered_map<std::string, TemplatePosition> &);
    void createDistributionDirectory();

    private:
    void openBamFile();

    public:
    Sample(Sample&&);
    Sample(std::string, SeqFileHandler &, ProgramOptions &, RegionSampler &, std::vector<complexVariant> &);
    void open();
    LibraryDistribution & getLibraryDistribution();
    void closeBamFile();
    int getFilterMargin();
    ContigInfo getContigInfo();
    std::unordered_map<std::string, int> getContigLengths();
    std::string getFileName();
    int getMaxReadLength();
    void close();
};

#endif
