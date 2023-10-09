#ifndef SAMPLEHEADER
#define SAMPLEHEADER

#include <fstream>
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
    ProgramOptions options;

    // local resources
    std::string filename;
    std::string sampleName;
    std::string distributionDirectory;
    std::string refFileName;
    BamFileHandler bamFile;
    SeqFileHandler referenceFile;

    std::vector<GenomicRegion> sampledRegions;
    std::unordered_map<std::string, int> chromosomeLengths;

    // results
    LibraryDistribution sampleDistribution;
    
    bool bamFileOpen;
    std::vector<std::string> regionStrings;
    int minMapQ;

    // private functions
    void calculateDefaultDistributions();
    void createSampleDistribution(std::unordered_map<std::string, TemplatePosition> &);
    void createDistributionDirectory();

    private:
    void openBamFile();

    public:
    Sample(Sample&&);
    Sample(std::string, std::string, ProgramOptions &, RegionSampler &);
    Sample(std::string);

    LibraryDistribution & getLibraryDistribution();
    void closeBamFile();
    int getFilterMargin();
    ContigInfo getContigInfo();
    std::unordered_map<std::string, int> getContigLengths();
    std::string getFileName();
    int getMaxReadLength();
    void close();
    void writeSampleProfile(std::ofstream &);
    void readSampleProfile(std::ifstream &);
    void printSampleProfile();
};

#endif
