#ifndef SAMPLEHEADER
#define SAMPLEHEADER

#include <fstream>
#include <vector>
#include <chrono>
#include <cstring>
#include <ostream>
#include <stdexcept>
#include <unordered_map>

#include "bamFileHandler.hpp"
#include "genotypeDistribution.hpp"
#include "genotypeResult.hpp"
#include "recordManager.hpp"
#include "regionSampler.hpp"
#include "seqFileHandler.hpp"
#include "variant.hpp"
#include "libraryDistribution.hpp"
#include "custom_types.hpp"
#include "filter.hpp"
#include "profileHandler.hpp"


class Sample
{
    std::string filename;
    std::string sampleName;
    std::string distributionDirectory;
    BamFileHandler bamFile;
    SeqFileHandler referenceFile;

    std::vector<GenomicRegion> sampledRegions;
    ContigInfo contigInfo;

    // results
    LibraryDistribution sampleDistribution;
    
    bool bamFileOpen;
    std::vector<std::string> regionStrings;
    int minMapQ;

    // private functions
    void calculateDefaultDistributions();

    private:
    void openBamFile();

    public:
    Sample(std::string, const std::vector<GenomicRegion> &);
    Sample(std::string);
    Sample();
    
    void loadPopDelProfile(PopDelProfileHandler &);

    LibraryDistribution & getLibraryDistribution();
    void closeBamFile();
    int getFilterMargin();
    const ContigInfo & getContigInfo();
    std::unordered_map<std::string, int32_t> getContigLengths();
    std::string getFileName();
    int getMaxReadLength();
    void close();
    void writeSampleProfile(std::string);
    void readSampleProfile(std::string);
    void printSampleProfile();
    std::string getSampleName();
};

#endif
