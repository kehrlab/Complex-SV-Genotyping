#ifndef LIBDISTHEADER
#define LIBDISTHEADER

#include "custom_types.hpp"
#include <eigen3/Eigen/Sparse>
#include "genomicRegion.hpp"
#include "readTemplate.hpp"
#include "seqFileHandler.hpp"
#include <vector>

class LibraryDistribution
{
    std::vector<std::vector<float>> distribution;
    std::vector<std::vector<float>> uniformDistribution;
    std::vector<std::vector<float>> correctionFactors;

    std::vector<int> insertSizes;
    std::vector<float> gcDistribution;
    std::vector<float> insertDistribution;

    float insertMean;
    float insertSD;
    bool gcCorrection;

    int calculateGCIndex(float);
    void smoothDistribution(std::vector<std::vector<float>> &, int);
    void scaleDistribution(std::vector<std::vector<float>> &);
    void createMarginalDistributions();
    void calculateCorrectionFactors();
    void calculateInsertStats();

    public:
    LibraryDistribution();
    LibraryDistribution(std::unordered_map<std::string, TemplatePosition> &);
    LibraryDistribution(std::unordered_map<std::string, TemplatePosition> &, std::vector<GenomicRegion> &, SeqFileHandler &);
    void writeDistribution(std::string filename);
    float getProbability(int, float);
    float getProbability(int);
    float getProbability(float);
    float getCorrectionFactor(int, float);
    float getInsertMean();
    float getInsertSD();
    bool gcCorrectionPossible();
    int drawInsertSize();
};

#endif