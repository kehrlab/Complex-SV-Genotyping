#include "libraryDistribution.hpp"
#include "variantProfile.hpp"

// 
// 
//

class GenotypeProfile
{
    std::vector<Eigen::SparseMatrix<float>> genotypeDistributions;
    std::vector<std::string> genotypeNames;
    std::unordered_map<std::string, int> genotypeIndices;
    std::unordered_map<std::string, int> groupIndices;
    int sMin;
    int sMax;

    public:
    GenotypeProfile(LibraryDistribution &, VariantProfile &);

    private:
    std::vector<Eigen::SparseMatrix<float>> calculateAlleleProfiles(LibraryDistribution &, VariantProfile &);
    void mixAlleleProfiles(std::vector<Eigen::SparseMatrix<float>>, VariantProfile &);
    void rescaleGenotypeProfiles();

    float getProbability(...);
    void writeDistribution();
};