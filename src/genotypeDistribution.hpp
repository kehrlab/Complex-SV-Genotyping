#ifndef GENOTYPEDISTRIBUTIONHEADER
#define GENOTYPEDISTRIBUTIONHEADER

#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <eigen3/Eigen/Sparse>
#include <cmath>
#include <stdexcept>

#include "insertSizeDistribution.hpp"

class GenotypeDistribution
{
    std::unordered_map<std::string, InsertSizeDistribution> distributions;
    std::unordered_map<std::string, float> interChromosomeCounts;
    std::unordered_map<std::string, float> distributionProbabilities;

    float normalizationFactor;
    std::unordered_set<std::string> contigs;

    float minProbability;
    int mode;
    void calculateMinProbability();
    void initDistributions();
    void initDistributionsDefault();
    void initDistributionsFine();
    void initDistributionsLocal();

    public:
    GenotypeDistribution();
    GenotypeDistribution(std::vector<std::string>, int);
    GenotypeDistribution(const Eigen::SparseMatrix<float> &, std::unordered_map<std::string, int> &, int, int);

    GenotypeDistribution & operator+=(GenotypeDistribution &);
    GenotypeDistribution & operator*=(float &);

    void addInsertSizeProbability(int, std::string, bool, bool, bool, std::string, std::string, std::string, std::vector<std::vector<std::string>>, float);
    void calculateNormalizationFactor();
    void scaleDistribution();

    static std::vector<std::string> determineDistributionKeys(std::string, bool, bool, int, std::string, std::string, std::string);
    float getProbability(int, std::string, bool, bool, bool, std::string, std::string, std::string, std::vector<std::vector<std::string>>, bool);

    void writeDistribution(std::string);
    void writeDistributionBinned(std::string);

    void setPossibleContigs(std::vector<std::string>);
    void setDistributionMode(int);

    std::string getContigIdentifier(std::string, std::string);
    void addInterChromosomeProbability(std::string, float);
    
    std::unordered_map<std::string, InsertSizeDistribution> & getDistributions();

    float getInterChromosomeProbability(std::string);
    std::unordered_map<std::string, float> getInterChromosomeProbabilites();
    std::unordered_set<std::string> getPossibleContigs();
    float getTotalInterChromosomeProbability();
    float getMinProbability();
    void setMinProbability(float);
    float calculateKLD(GenotypeDistribution &);
    int getMinInsertSize();
    int getMaxInsertSize();
    void smoothDistribution();
    float getNormalizationFactor();
};

GenotypeDistribution operator+(GenotypeDistribution, GenotypeDistribution);
GenotypeDistribution operator*(float, GenotypeDistribution);
GenotypeDistribution operator*(GenotypeDistribution, float);

#endif
