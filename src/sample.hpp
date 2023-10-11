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
#include "libraryDistribution.hpp"
#include <chrono>

class Sample
{
    std::string filename;
    std::string sampleName;
    std::string distributionDirectory;
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
    // void createSampleDistribution(std::unordered_map<std::string, TemplatePosition> &);
    // void createDistributionDirectory();

    private:
    void openBamFile();

    public:
    Sample & operator=(Sample);
    // Sample(Sample&&);
    Sample(std::string, const std::vector<GenomicRegion> &);
    Sample(std::string);
    

    LibraryDistribution & getLibraryDistribution();
    void closeBamFile();
    int getFilterMargin();
    ContigInfo getContigInfo();
    std::unordered_map<std::string, int> getContigLengths();
    std::string getFileName();
    int getMaxReadLength();
    void close();
    void writeSampleProfile(std::string);
    void readSampleProfile(std::string);
    void printSampleProfile();
    std::string getSampleName();
};

#endif
