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
    float minProbability;
    

    public:
    GenotypeDistribution();
    GenotypeDistribution(const Eigen::SparseMatrix<float> &, std::unordered_map<std::string, int> &, int, int);

    GenotypeDistribution & operator+=(GenotypeDistribution &);
    GenotypeDistribution & operator*=(float &);

    void addInsertSizeProbability(int, std::string, std::string, std::string, std::string, float);
    void addInterChromosomeProbability(std::string, float);

    void calculateNormalizationFactor();
    void scaleDistribution();

    float getProbability(int, std::string, std::string, std::string, std::string, bool &);
    float getInterChromosomeProbability(std::string);

    void writeDistribution(std::string);
    void writeDistributionBinned(std::string);
    
    std::unordered_map<std::string, InsertSizeDistribution> & getDistributions();
    std::unordered_map<std::string, float> getInterChromosomeProbabilites();
    float getTotalInterChromosomeProbability();

    void calculateMinProbability();
    float getMinProbability();
    void setMinProbability(float);
    float calculateKLD(GenotypeDistribution &);
    int getMinInsertSize();
    int getMaxInsertSize();
    void smoothDistribution();
    float getNormalizationFactor();

    static std::string determineGroup(std::string orientation, std::string junctionString, std::string breakpointString, std::string chromosomeString)
    {
        std::string group;
        if (junctionString != "")
            group = "split_" + junctionString;
        else if (breakpointString != "")
            group = "spanning_" + breakpointString;
        else if (orientation != "")
            group = orientation;
        else 
        {
            std::cerr << "Could not assign a group to read pair." << std::endl;
            group = "";
        }
        if (group == "FR")
            group = "RF";
        return group;
    }
};

GenotypeDistribution operator+(GenotypeDistribution, GenotypeDistribution);
GenotypeDistribution operator*(float, GenotypeDistribution);
GenotypeDistribution operator*(GenotypeDistribution, float);

#endif
