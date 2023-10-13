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

    /// 
    std::vector<float> genotypeLikelihoods;
    std::vector<float> genotypeLikelihoodSDs;
    std::vector<std::string> genotypeNames;
    std::vector<int> mappingQualities;
    std::unordered_map<std::string, std::vector<int>> readPairGroups;
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
    int mode;
    std::string calledGenotype;
    int minQuality;
    int maxQuality;
    float meanQuality;
    std::unordered_map<std::string, float> groupQualities;

    void getLikelihoodsFromBootstrapping();
    void determineQualityStats();
    void addPairGroups(std::vector<std::string>);

    public:
    GenotypeResult();
    GenotypeResult(std::string, bool);
    GenotypeResult(std::string, std::string, std::vector<std::string>, bool);
    
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
    void storeEvidence(int, std::string, bool, bool, bool, std::string, std::string, std::string, std::vector<std::vector<std::string>>, std::vector<int>);
    void writeEvidence(std::string);
    void setPossibleContigs(std::vector<std::string>);
    void setDistributionMode(int);
    void printAllLikelihoods();
    void addTemplateProbabilities(std::vector<std::string>, std::vector<float>, float);
    void addProbability(std::string, float);
    void calculateLikelihoods();
    void bootstrapLikelihoods();
    float mean(std::vector<float>);
    float sd(std::vector<float>, float);
    void writeBootstrapData(std::string);
    void bootstrapQuality();
    void clearData();
    void calculateReadStats();
    std::unordered_map<std::string, float> getReadStats();
    std::string getSampleName();
};

#endif