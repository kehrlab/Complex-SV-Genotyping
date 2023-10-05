#ifndef VARPROFILEHEADER
#define VARPROFILEHEADER

#include "filter.hpp"
#include "libraryDistribution.hpp"
#include "options.hpp"
#include "variant.hpp"
#include <eigen3/Eigen/Sparse>
#include <unordered_map>
#include "genomicRegion.hpp"
#include "genotypeDistribution.hpp"
#include "custom_types.hpp"

class VariantProfile
{
    complexVariant variant;
    ProgramOptions options;
    std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> chromosomeStructures; // need one of these for every variant allele!

    std::unordered_map<std::string, int> variantAlleleNames;
    std::unordered_map<std::string, int> variantGroups;

    Eigen::SparseMatrix<float, Eigen::RowMajor> referenceMask;
    std::vector<std::vector<Eigen::SparseMatrix<float, Eigen::RowMajor>>> variantMask;

       
    int sMin;
    int sMax;

    int sMinMapped;
    int sMaxMapped;

    int filterMargin;
    int overlap;
    int readLength;

    ReadPairFilter filter;
    

    public:
    VariantProfile();
    VariantProfile(complexVariant, int, int, int, int, int, const std::unordered_map<std::string, int> &, ProgramOptions &);

    void calculateAlleleMasks();

    const Eigen::SparseMatrix<float, Eigen::RowMajor> & getVariantMask(int s);
    const Eigen::SparseMatrix<float, Eigen::RowMajor> & getVariantMask(const std::string &, int s);
    const Eigen::SparseMatrix<float, Eigen::RowMajor> & getReferenceMask();

    std::unordered_map<std::string, GenotypeDistribution> calculateGenotypeDistributions(LibraryDistribution &, float);
    ReadPairFilter & getFilter();
    complexVariant & getVariant();
    std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> & getChromosomeStructures();

    private:
    void determinePossibleGroups();
    std::vector<std::unordered_set<int>> createIndexCombinations(std::vector<int> &, std::vector<int> &);
    std::vector<std::unordered_set<int>> createIndexCombinations(std::vector<int> &, std::vector<int> &, std::unordered_map<int, std::unordered_set<int>> &);
    std::unordered_map<int, std::unordered_set<int>> getJunctionAmbiguities(VariantRegions &);

    inline void createIndexString(std::string &, const std::unordered_set<int> &);
    inline void createIndexString(std::string &, const std::unordered_set<std::string> &);
    std::vector<GenomicRegion> getRegionsFromIndices(JunctionRegion &, int, int);
    void insertBreakpoint(Breakpoint &, JunctionRegion &);
    VariantRegions createVariantRegions(Allele &);
    void createChromosomeStructures(Allele &);

    void determineVariantGroups(VariantRegions &);
    void findPairAttributes(std::unordered_set<std::string> & groups, VariantRegions & variantRegions);
    void determineSpanningGroups(std::unordered_set<std::string> & groups, VariantRegions & variantRegions);
    void determineSplitGroups(std::unordered_set<std::string> & groups, VariantRegions & variantRegions);

    void initMasks();
    void initReferenceMask();
    void initVariantMask();
    inline void addSimulatedTemplateToMask(int &, VariantMap &, int, int, Allele &);
    inline void addValueToMask(Allele & allele, int sOld, int sNew, std::string & orientation, std::string & jString, std::string & bpString, std::string & chromosomes);
};

#endif