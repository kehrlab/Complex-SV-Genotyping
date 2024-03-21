#ifndef INSERTSIZEDISTHEADER
#define INSERTSIZEDISTHEADER

#include <random>
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <eigen3/Eigen/Sparse>

#include "breakpoint.hpp"
#include "junction.hpp"


class InsertSizeDistribution
{
    double insertMean;
    double insertSD;

    int64_t minInsertSize;
    int64_t maxInsertSize;

    Eigen::SparseVector<float> distribution;

    float minValue;

    public: 
    InsertSizeDistribution();
    InsertSizeDistribution(std::vector<int64_t>);
    InsertSizeDistribution(int64_t, int64_t, Eigen::SparseVector<float>);

    InsertSizeDistribution & operator+=(InsertSizeDistribution & rhs);
    InsertSizeDistribution & operator*=(float & factor);

    std::vector<int> variantToReferenceMap(std::vector<Junction>, int);
    std::vector<int> getJunctionIndices(std::vector<int>);

    void addInsertSizeProbability(int64_t, float);
    void setMinInsertSize(int);
    void setMaxInsertSize(int);
    void inferMinAndMaxValues(std::vector<int64_t> const &);
    void generateDistributionFromInsertSizes(std::vector<int64_t> const &);
    void calculateDistributionMean(std::vector<int64_t> const &);
    void calculateDistributionSD(std::vector<int64_t> const &);
    void scaleDistribution();
    void smoothDistribution(int);
    void writeDistribution(std::string);
    double getInsertMean() const;
    double getInsertSD() const;
    float getInsertSizeProbability(int64_t);
    int64_t getMinInsertSize();
    int64_t getMaxInsertSize();
    float getIntegral();
    Eigen::SparseVector<float> & getDistributionVector();
    void findMinProbability();
    float getMinProbability();
};

InsertSizeDistribution operator+(InsertSizeDistribution, InsertSizeDistribution);
InsertSizeDistribution operator*(float, InsertSizeDistribution);
InsertSizeDistribution operator*(InsertSizeDistribution, float);

#endif