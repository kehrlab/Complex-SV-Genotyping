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
    std::unordered_map<std::string, float> distributionProbabilities;

    float normalizationFactor;
    float minProbability;
    

    public:
    GenotypeDistribution();
    GenotypeDistribution(const Eigen::SparseMatrix<float, Eigen::RowMajor> &, std::unordered_map<std::string, int> &, int64_t, int64_t);

    GenotypeDistribution & operator+=(GenotypeDistribution &);
    GenotypeDistribution & operator*=(float &);

    void addInsertSizeProbability(int64_t, std::string, std::string, std::string, std::string, float);
    void calculateNormalizationFactor();
    void scaleDistribution();
    float getProbability(int64_t, std::string, std::string, std::string, std::string, bool &);
    void writeDistribution(std::string);
    void writeDistributionBinned(std::string);

    std::unordered_map<std::string, InsertSizeDistribution> & getDistributions();

    void calculateMinProbability();
    float getMinProbability();
    void setMinProbability(float);
    float calculateKLD(GenotypeDistribution &);
    int64_t getMinInsertSize();
    int64_t getMaxInsertSize();
    void smoothDistribution();
    float getNormalizationFactor();

    static std::string determineGroup(std::string orientation, std::string junctionString, std::string breakpointString, std::string bridgeString)
    {
        std::string group;

        if (junctionString != "")
            group = "split_" + junctionString;
        else if (breakpointString != "")
            group = "spanning_" + breakpointString;
        // else if (bridgeString != "")
        //     group = "bridging_" + bridgeString;
        else if (orientation != "")
            group = orientation;
        else 
        {
            std::cerr << "Could not assign a group to read pair." << std::endl;
            std::cerr << "Orientation: " << orientation << std::endl;
            std::cerr << "Junctions: " << junctionString << std::endl;
            std::cerr << "Breakpoints: " << breakpointString << std::endl;
            std::cerr << "Bridged breakpoints: " << bridgeString << std::endl;
            std::cerr << std::endl;
            group = "";
        }
        if (group == "RF")
            group = "FR";
        return group;
    }
};

GenotypeDistribution operator+(GenotypeDistribution, GenotypeDistribution);
GenotypeDistribution operator*(float, GenotypeDistribution);
GenotypeDistribution operator*(GenotypeDistribution, float);

#endif
