#ifndef INSERTSIZEDISTHEADER
#define INSERTSIZEDISTHEADER

#include <vector>
#include <cmath>
#include <iostream>
#include "breakpoint.hpp"
#include "junction.hpp"
#include <random>

class InsertSizeDistribution
{
    double insertMean;
    double insertSD;

    int minInsertSize;
    int maxInsertSize;

    std::vector<float> distribution;

    float minValue;

    public: 
    InsertSizeDistribution();
    InsertSizeDistribution(std::vector<int>);
    InsertSizeDistribution(int, int, std::vector<float>);

    InsertSizeDistribution & operator+=(InsertSizeDistribution & rhs);
    InsertSizeDistribution & operator*=(float & factor);

    std::vector<int> variantToReferenceMap(std::vector<Junction>, int);
    std::vector<int> getJunctionIndices(std::vector<int>);

    void addInsertSizeProbability(int, float);
    void setMinInsertSize(int);
    void setMaxInsertSize(int);
    void inferMinAndMaxValues(std::vector<int> const &);
    void generateDistributionFromInsertSizes(std::vector<int> const &);
    int getNumberOfInsertSizes(std::vector<int> const &, int);
    void calculateDistributionMean(std::vector<int> const &);
    void calculateDistributionSD(std::vector<int> const &);
    void scaleDistribution();
    void smoothDistribution(int);
    void writeDistribution(std::string);
    double getInsertMean() const;
    double getInsertSD() const;
    float getInsertSizeProbability(int);
    int getMinInsertSize();
    int getMaxInsertSize();
    float getIntegral();
    void findMinProbability();
    float getMinProbability();
};

InsertSizeDistribution operator+(InsertSizeDistribution, InsertSizeDistribution);
InsertSizeDistribution operator*(float, InsertSizeDistribution);
InsertSizeDistribution operator*(InsertSizeDistribution, float);

#endif