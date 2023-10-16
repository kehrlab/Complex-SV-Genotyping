#ifndef VARIANTHEADER
#define VARIANTHEADER

#include <iterator>
#include <iostream>
#include <vector>
#include <unordered_set>

#include "genomicRegion.hpp"
#include "insertSizeDistribution.hpp"
#include "libraryDistribution.hpp"
#include "allele.hpp"
#include "variantParser.hpp"
#include "breakpoint.hpp"
#include "seqFileHandler.hpp"


class complexVariant {
    std::vector<Allele> variantAlleles;
    std::string variantName;
    std::vector<Breakpoint> allBreakpoints;
    std::vector<Junction> allJunctions;

    std::vector<GenomicRegion> variantRegions;

    int nAlleles;
    int filterMargin;
    std::string variantFileName;

    int calculateSearchDistance(double, double);
    void poolBreakpoints();
    void poolJunctions();

    std::vector<GenomicRegion> createRegionsFromBreakpoints(int);

    public:
    complexVariant();
    complexVariant(std::string, std::vector<std::string>, variantData);
    complexVariant(std::string, std::vector<std::string>, variantData, std::string);
    complexVariant(std::vector<Junction>);

    void print();
    std::vector<GenomicRegion> calculateAssociatedRegions(LibraryDistribution &);
    std::string getName();
    std::vector<Junction> & getAllJunctions();
    std::vector<Breakpoint> & getAllBreakpoints();
    std::vector<Allele> & getAlleles();
    std::unordered_set<std::string> getContigNames();
    void calculateVariantRegions();
    std::vector<GenomicRegion> & getVariantRegions();
    void createAlleleMaps(int, std::unordered_map<std::string, int>);
    std::string getVariantFileName();
    void setFilterMargin(int);
};

#endif