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

    // ProgramOptions options;
    std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> chromosomeStructures;

    std::unordered_map<std::string, int> variantAlleleNames;
    std::unordered_map<std::string, int> variantGroups;

    Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t> referenceMask;
    std::vector<std::vector<Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t>>> variantMask;
       
    int filterMargin;
    int overlap;
    int readLength;

    int sMin;
    int sMax;

    ContigInfo cInfo;

    int64_t sMinMapped;
    int64_t sMaxMapped;
    

    ReadPairFilter filter;
    

    public:
    VariantProfile();
    VariantProfile(std::string);
    VariantProfile(int, int, int, int, int, const ContigInfo &);
    VariantProfile(complexVariant, int, int, int, int, int, const ContigInfo &);

    void calculateAlleleMasks();

    const Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t> & getVariantMask(int s);
    const Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t> & getVariantMask(const std::string &, int s);
    const Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t> & getReferenceMask();

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
    void createVariantChromosomeStructures();
    std::string getName() const;
    const ContigInfo & getContigInfo() const;
    void clearProfile();
    void overwrite(complexVariant);

    private:
    void determinePossibleGroups();

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
    bool isPossible(std::vector<int> positions);
    std::vector<std::vector<int>> getSubsets(std::vector<int>, std::vector<int>, uint32_t);

    void initMasks();
    void resizeMasks();
    void initReferenceMask();
    void initVariantMask();
    inline void addSimulatedTemplateToMask(int &, VariantMap &, int, int64_t, Allele &);
    inline void addValueToMask(Allele & allele, int64_t sOld, int64_t sNew, std::string & orientation, std::string & jString, std::string & bpString, std::string & chromosomes);
    void addGroupToMasks(std::string);

    inline void readString(std::string &, std::ifstream &);
    inline void writeString(std::ofstream &, std::string);
    void writeJunctions(std::ofstream &, std::string, std::vector<Junction> &);
    std::vector<Junction> readJunctions(std::ifstream &);
};

#endif
