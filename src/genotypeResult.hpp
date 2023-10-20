#ifndef GENOTYPERESULTHEADER
#define GENOTYPERESULTHEADER

#include <iostream>
#include <random>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "genotypeDistribution.hpp"

class GenotypeResult
{
    std::string filename;
    std::string sampleName;
    bool useQualities;
    int observedReads;
    int outlierCount;

    /// 
    std::vector<float> genotypeLikelihoods;
    std::vector<float> genotypeLikelihoodSDs;
    std::vector<std::string> genotypeNames;
    std::vector<int> mappingQualities;

    std::vector<std::vector<float>> templateProbabilities;
    std::vector<float> templateWeights;
    std::vector<std::vector<float>> bootstrappedLikelihoods;
    std::vector<float> bootstrappedQualities;

    GenotypeDistribution sampleDistribution;

    // result
    std::string outputString;
    float quality;
    float lowerBoundQuality;
    float upperBoundQuality;
    float callCertainty;

    std::string calledGenotype;
    int minQuality;
    int maxQuality;
    float meanQuality;

    void determineQualityStats();

    public:
    GenotypeResult();
    GenotypeResult(std::string, bool);
    GenotypeResult(std::string, std::string, bool);
    
    void callGenotype();
    void scaleToPhred();
    void initLikelihoods(std::vector<std::string>);
    float getLikelihood(std::string);
    std::string getCalledGenotype();
    void setFilename(std::string);
    std::string getFilename();
    float getQuality();
    void createOutputString();
    std::string getOutputString();
    void storeEvidence(int, std::string, std::string, std::string, std::string, std::vector<int>);
    void writeEvidence(std::string);
    void printAllLikelihoods();
    void addTemplateProbabilities(std::vector<std::string>, std::vector<float>, float);
    void addProbability(std::string, float);
    void calculateLikelihoods();
    void bootstrapLikelihoods();
    float mean(std::vector<float>);
    float sd(std::vector<float>, float);
    void writeBootstrapData(std::string);
    void bootstrapQuality(int, int);
    void clearData();
    void addOutlier(bool);
    std::string getSampleName();
};

#endif