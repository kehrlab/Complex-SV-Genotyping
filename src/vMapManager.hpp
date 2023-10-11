#ifndef VMAPMANAGERHEADER
#define VMAPMANAGERHEADER

#include <unordered_map>
#include <vector>
#include <string>

#include "breakpoint.hpp"
#include "filter.hpp"
#include "junction.hpp"
#include "variantMap.hpp"
#include "genomicRegion.hpp"


class VariantMapManager
{
    VariantMap variantMap;
    std::vector<VariantMap> refMaps;

    GenomicRegion refRegionLeft;
    GenomicRegion refRegionRight;

    int filterMargin;
    std::unordered_map<std::string, int> contigLengths;
    std::string chromosome;

    bool variantMapPresent;

    void setFilterMargin(int);
    void setContigLengths(std::unordered_map<std::string, int>);
    void setChromosomeName(std::string);

    void mergeOverlappingMaps();
    bool isOverlap(VariantMap &, VariantMap &);

    void fillMap(VariantMap &, Breakpoint &);
    void fillMap(VariantMap &, std::vector<Junction> &);
    void fillPositionsLeft(VariantMap &, std::vector<Junction> &);
    void fillPositionsBetweenAllJunctions(VariantMap &, std::vector<Junction> &);
    void fillPositionsBetweenJunctions(VariantMap &, Junction &, Junction &, int);
    void fillPositionsRight(VariantMap &, std::vector<Junction> &);

    void createReferenceRegions();

    void extendMapRight(VariantMap &, GenomicRegion &);
    void extendMapLeft(VariantMap &, GenomicRegion &);

    public:
    VariantMapManager();
    VariantMapManager(std::string, int, std::unordered_map<std::string, int>);
    VariantMapManager(std::string, std::vector<Junction>, std::vector<Breakpoint>, int, std::unordered_map<std::string, int>);

    void createMapFromJunctions(std::vector<Junction>);
    void createBreakpointMaps(std::vector<Breakpoint>);
    void createMapFromBreakpoint(Breakpoint);

    std::vector<VariantMap> getMaps();
    void print();
    seqan::String<seqan::Dna5String> getChromosomeSequences(SeqFileHandler &);
};

#endif