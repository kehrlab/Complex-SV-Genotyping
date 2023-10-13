#ifndef VARPROFILEHEADER
#define VARPROFILEHEADER

#include <eigen3/Eigen/Sparse>
#include <unordered_map>
#include <fstream>
#include <ios>
#include <stdexcept>
#include <unordered_set>

#include "filter.hpp"
#include "libraryDistribution.hpp"
#include "variant.hpp"
#include "genomicRegion.hpp"
#include "genotypeDistribution.hpp"
#include "custom_types.hpp"
#include "readTemplate.hpp"
#include "variantParser.hpp"


class VariantProfile
{
    complexVariant variant;
    std::string name;
    bool variantPresent;
    // ProgramOptions options;
    std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> chromosomeStructures;

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
    VariantProfile(std::string);
    VariantProfile(complexVariant, int, int, int, int, int, const std::unordered_map<std::string, int> &);

    void calculateAlleleMasks();

    const Eigen::SparseMatrix<float, Eigen::RowMajor> & getVariantMask(int s);
    const Eigen::SparseMatrix<float, Eigen::RowMajor> & getVariantMask(const std::string &, int s);
    const Eigen::SparseMatrix<float, Eigen::RowMajor> & getReferenceMask();

    void calculateGenotypeDistributions(std::unordered_map<std::string, GenotypeDistribution> &, LibraryDistribution &, float);
    ReadPairFilter & getFilter();
    complexVariant & getVariant();
    std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> & getChromosomeStructures();

    void writeProfile(std::string);
    void readProfile(std::string);
    int getMinInsert();
    int getMaxInsert();
    int getMargin();
    int getReadLength();
    bool variantStructureIsPresent();
    bool loadVariantStructure(std::string filename, std::string variantName);
    std::string getName();

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