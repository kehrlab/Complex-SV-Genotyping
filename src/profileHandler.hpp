#ifndef PROFILEHANDLER
#define PROFILEHANDLER

#include <string>
#include <unordered_map>
#include "genomicRegion.hpp"
#include "libraryDistribution.hpp"
#include "readTemplate.hpp"

#include "popdel/popdel_snippets.hpp"


class PopDelProfileHandler{

    private:
        std::string profilePath;
        std::string sampleName;
        std::ifstream infile;

        unsigned sampleNumRegions;
        unsigned indexRegionSize;
        unsigned numReadGroups;
        seqan::String<seqan::CharString> sampleContigNames;
        seqan::String<int32_t> sampleContigLengths;

        std::vector<std::string> readGroups;
        std::vector<Histogram> histograms; // need to integrate PopDel code

        LibraryDistribution sampleDistribution;

    public:
        PopDelProfileHandler();
        PopDelProfileHandler(std::string);

        LibraryDistribution & getLibraryDistribution();
        std::string getSampleName();
        std::string getFileName();

        std::vector<std::string> getContigNames();
        std::vector<int> getContigLengths();

        std::unordered_map<std::string, std::vector<BamRecord>> get_read_pairs(std::vector<GenomicRegion>);
        inline void createReadPair(std::vector<BamRecord> &, const ReadPair &, seqan::CharString, int);
        inline void createReadPair(std::vector<BamRecord> &, const TranslocationWindowEntry &, seqan::CharString, int);

        void open(std::string);
        void close();

    private:
        inline void readHistograms();
        inline void mergeRGHistograms();
};

#endif