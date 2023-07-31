#include "genotypeResult.hpp"
#include "genotypeDistribution.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>


GenotypeResult::GenotypeResult()
{
    this->filename = "";
    this->observedReads = 0;
    this->callCertainty = 0;
    this->lowerBoundQuality = 0;
    this->upperBoundQuality = 0;
    this->minQuality = 0;
    this->maxQuality = 0;
    this->meanQuality = 0.0;
    this->useQualities = false;
}

GenotypeResult::GenotypeResult(std::string filename, bool useQualities)
{
    this->filename = filename;
    this->observedReads = 0;
    this->callCertainty = 0;
    this->lowerBoundQuality = 0;
    this->upperBoundQuality = 0;
    this->minQuality = 60;
    this->maxQuality = 0;
    this->meanQuality = 0.0;
    this->useQualities = useQualities;
}

GenotypeResult::GenotypeResult(std::string filename, std::vector<std::string> contigNames, bool useQualities)
{
    this->filename = filename;
    this->observedReads = 0;
    setPossibleContigs(contigNames);
    this->callCertainty = 0;
    this->lowerBoundQuality = 0;
    this->upperBoundQuality = 0;
    this->minQuality = 60;
    this->maxQuality = 0;
    this->meanQuality = 0.0;
    this->useQualities = useQualities;
}

GenotypeResult::GenotypeResult(std::string filename, std::vector<std::string> contigNames, ProgramOptions & options)
{
    this->filename = filename;
    this->observedReads = 0;
    this->callCertainty = 0;
    this->lowerBoundQuality = 0;
    this->upperBoundQuality = 0;
    this->minQuality = 60;
    this->maxQuality = 0;
    this->meanQuality = 0.0;
    this->useQualities = options.isOptionUseQualities();
    setPossibleContigs(contigNames);
    setDistributionMode(options.getDistributionMode());
}

void GenotypeResult::initLikelihoods(std::vector<std::string> genotypeNames)
{
    this->genotypeNames = genotypeNames;
    for (auto gtName : genotypeNames) {
        std::vector<float> emptyVec;
        this->templateProbabilities.push_back(emptyVec);
    }
}

void GenotypeResult::callGenotype()
{
    bootstrapLikelihoods();
    getLikelihoodsFromBootstrapping();
    bootstrapQuality();
    

    if (this->genotypeLikelihoods.size() < 3)
        throw std::runtime_error("ERROR: There must be at least 3 different genotypes.");

    float lMin = 0;
    int minIndex;

    std::vector<float>::iterator result = std::min_element(this->genotypeLikelihoods.begin(), this->genotypeLikelihoods.end());
    lMin = *result;
    minIndex = std::distance(this->genotypeLikelihoods.begin(), result);

    float lSecond = *(std::max_element(this->genotypeLikelihoods.begin(), this->genotypeLikelihoods.end()));
    for (int i = 0; i < this->genotypeLikelihoods.size(); ++i)
    {
        if (i == minIndex)
            continue;
        lSecond = std::min(lSecond, this->genotypeLikelihoods[i]);
    }

    this->calledGenotype = this->genotypeNames[minIndex];
    this->quality = lSecond - lMin;

    calculateReadStats();
    createOutputString();
}

std::string GenotypeResult::getCalledGenotype()
{
    return this->calledGenotype;
}

std::string GenotypeResult::getFilename()
{
    return this->filename;
}

void GenotypeResult::setFilename(std::string filename)
{
    this->filename = filename;
}

void GenotypeResult::addTemplateProbabilities(std::vector<std::string> genotypeNames, std::vector<float> probabilities, float weight)
{
    this->templateWeights.push_back(weight);
    for (int i = 0; i < genotypeNames.size(); ++i)
        addProbability(genotypeNames[i], probabilities[i]);
    return;
}

void GenotypeResult::addProbability(std::string genotypeName, float probability)
{
    int gtIdx = -1;
    for (int i = 0; i < this->genotypeNames.size(); ++i)
    {
        if (this->genotypeNames[i] == genotypeName)
        {
            gtIdx = i;
            break;
        }
    }
    if (gtIdx != -1)
        this->templateProbabilities[gtIdx].push_back((-10) * std::log10(probability));
    else
        throw std::runtime_error("Did not find Genotype to adjust!");
    return;
}

void GenotypeResult::calculateLikelihoods()
{
    for (int i = 0; i < this->genotypeNames.size(); ++i)
    {
        float L = 0;
        for (int j = 0; j < this->templateProbabilities[i].size(); ++j)
            L += this->templateProbabilities[i][j];
        this->genotypeLikelihoods.push_back(L);
    }
}

void GenotypeResult::bootstrapLikelihoods()
{
    std::random_device rd;
    std::mt19937 gen(rd());


    int n = this->templateProbabilities[0].size();
    int sampleSize = n;
    int nSamples = 10000;
    int idx = -1;

    if (!this->useQualities)
        for (float & w : this->templateWeights)
            w = 1.0;
        
    std::discrete_distribution<> distrib(this->templateWeights.begin(), this->templateWeights.end());

    // init likelihood vector
    for (int i = 0; i < this->genotypeNames.size(); ++i)
    {
        std::vector<float> temp;
        this->bootstrappedLikelihoods.push_back(temp);   
    }
    
    float L = 0.0;
    for (unsigned j = 0; j < nSamples; ++j)
    {
        // create bootstrap sample
        std::vector<int> indices(sampleSize, 0);
        for (unsigned k = 0; k < sampleSize; ++k)
            indices[k] = distrib(gen);

        for (int i = 0; i < this->genotypeNames.size(); ++i)
        {
            L = 0.0;
            for (unsigned k = 0; k < sampleSize; ++k) 
                L += this->templateProbabilities[i][indices[k]];
            this->bootstrappedLikelihoods[i].push_back(L);
        }
    }
}

void GenotypeResult::getLikelihoodsFromBootstrapping()
{
    this->genotypeLikelihoods.erase(this->genotypeLikelihoods.begin(), this->genotypeLikelihoods.end());
    for (int i = 0; i < this->genotypeNames.size(); ++i)
    {
        float gtMean = mean(this->bootstrappedLikelihoods[i]);
        float gtSD = sd(this->bootstrappedLikelihoods[i], gtMean);
        this->genotypeLikelihoods.push_back(gtMean);
        this->genotypeLikelihoodSDs.push_back(gtSD);
    }
}

void GenotypeResult::bootstrapQuality()
{
    std::vector<float> likelihoods = this->genotypeLikelihoods;
    std::sort(likelihoods.begin(), likelihoods.end());
    int minIdx = std::find(this->genotypeLikelihoods.begin(), this->genotypeLikelihoods.end(), likelihoods[0]) - this->genotypeLikelihoods.begin();
    int secondIdx = std::find(this->genotypeLikelihoods.begin(), this->genotypeLikelihoods.end(), likelihoods[1]) - this->genotypeLikelihoods.begin();

    int nSameCall = 0;
    float q = 0;
    for (int i = 0; i < this->bootstrappedLikelihoods[minIdx].size(); ++i)
    {
        q = this->bootstrappedLikelihoods[secondIdx][i] - this->bootstrappedLikelihoods[minIdx][i];
        this->bootstrappedQualities.push_back(q);
        if (q > 0)
            ++nSameCall;
    }
    this->callCertainty = (float) nSameCall / (float) this->bootstrappedLikelihoods[minIdx].size();

    std::vector<float> bQs = this->bootstrappedQualities;
    std::sort(bQs.begin(), bQs.end());
    int minPercentileIdx = (int) (0.025 * bQs.size()) - 1;
    int maxPercentileIdx = (int) (0.975 * bQs.size()) - 1;
    this->lowerBoundQuality = bQs[minPercentileIdx];
    this->upperBoundQuality = bQs[maxPercentileIdx];
}

void GenotypeResult::writeBootstrapData(std::string prefix)
{
    std::string bsFilename = prefix + "_bootstrap_data.txt";
    std::ofstream f(bsFilename);
    if (f.is_open())
    {
        for (auto gt : this->genotypeNames)
            f << gt << "\t";
        f << "Quality" << std::endl;

        for (int i = 0; i < this->bootstrappedLikelihoods[0].size(); ++i)
        {
            std::vector<float> v;
            float minL = this->bootstrappedLikelihoods[0][i];
            for (int j = 0; j < this->bootstrappedLikelihoods.size(); ++j) {
                float L = this->bootstrappedLikelihoods[j][i];
                v.push_back(L);
                minL =  std::min(minL, L);
            }
            for (int j = 0; j < v.size(); ++j)
                f << (v[j] - minL) << "\t";
            f << this->bootstrappedQualities[i] << std::endl;
        }
    }
    f.close();
}

float GenotypeResult::mean(std::vector<float> v)
{
    double m = 0;
    for (float x : v)
        m += x;
    m /= v.size();
    return (float) m;
}

float GenotypeResult::sd(std::vector<float> v, float m)
{
    double s = 0;
    for (float x : v)
        s += (x-m)*(x-m);
    s /= v.size();
    s = std::sqrt(s);
    return (float) s;
}

float GenotypeResult::getLikelihood(std::string genotype)
{
    int index = -1;
    for (int i = 0; i < this->genotypeNames.size(); ++i)
    {
        if (this->genotypeNames[i] == genotype)
        {
            index = i;
            break;
        }
    }
    if (index == -1)
        throw std::runtime_error("Did not find Genotype to adjust!");
    return this->genotypeLikelihoods[index];
}

void GenotypeResult::printAllLikelihoods()
{
    std::cout << std::endl;
    for (int i = 0; i < this->genotypeNames.size(); ++i)
        std::cout << this->genotypeNames[i] << "\t";
    std::cout << std::endl;
    for (int i = 0; i < this->genotypeLikelihoods.size(); ++i)
        std::cout << this->genotypeLikelihoods[i] << "\t";
    std::cout << std::endl << std::endl;
}

void GenotypeResult::createOutputString()
{
    determineQualityStats();
    std::string outputString = "";
    outputString.append(this->filename);
    outputString.append("\t");
    outputString.append(this->calledGenotype);
    outputString.append("\t");
    outputString.append(std::to_string(this->quality));
    outputString.append("\t");
    outputString.append(std::to_string(this->lowerBoundQuality));
    outputString.append("\t");
    outputString.append(std::to_string(this->upperBoundQuality));
    outputString.append("\t");
    outputString.append(std::to_string(this->callCertainty));
    outputString.append("\t");
    outputString.append(std::to_string(this->observedReads));
    outputString.append("\t");
    outputString.append(std::to_string(this->meanQuality));
    outputString.append("\t");
    outputString.append(std::to_string(this->minQuality));
    outputString.append("\t");
    outputString.append(std::to_string(this->maxQuality));
    outputString.append("\t");
    if (this->meanQuality >= 55 && this->callCertainty >= 0.9999)
	    outputString.append("TRUE");
    else
	    outputString.append("FALSE");
    outputString.append("\n");
    this->outputString = outputString;
}

std::string GenotypeResult::getOutputString()
{
    return this->outputString;
}

void GenotypeResult::determineQualityStats()
{
    for (auto mapQ : this->mappingQualities)
    {
        if (mapQ < this->minQuality)
            this->minQuality = mapQ;
        if (mapQ > this->maxQuality)
            this->maxQuality = mapQ;
        this->meanQuality += (float) mapQ;
    }
    this->meanQuality /= this->mappingQualities.size();
}

void GenotypeResult::storeEvidence(int insertSize, std::string orientation, bool split, bool spanning, bool interChromosome, std::string regionString, std::string junctionString, std::string breakpointString, std::vector<std::vector<std::string>> rNamePairs, std::vector<int> mapQs)
{
    ++ this->observedReads;
    for (int mapQ : mapQs)
        this->mappingQualities.push_back(mapQ);
    std::vector<std::string> keys = GenotypeDistribution::determineDistributionKeys(orientation, split, spanning, mode, regionString, junctionString, breakpointString);
    addPairGroups(keys);
    this->sampleDistribution.addInsertSizeProbability(insertSize, orientation, split, spanning, interChromosome, regionString, junctionString, breakpointString, rNamePairs, 1.0);
}

void GenotypeResult::addPairGroups(std::vector<std::string> keys)
{
    for (auto key : keys) {
        if (this->readPairGroups.find(key) != this->readPairGroups.end())
            this->readPairGroups[key].push_back(this->observedReads - 1);
        else
            this->readPairGroups[key] = std::vector<int>(1, this->observedReads - 1);
    }
}

void GenotypeResult::writeEvidence(std::string prefix)
{
    this->sampleDistribution.scaleDistribution();
    this->sampleDistribution.writeDistributionBinned(prefix);
}

void GenotypeResult::setPossibleContigs(std::vector<std::string> contigNames)
{
    this->sampleDistribution.setPossibleContigs(contigNames);
}

void GenotypeResult::setDistributionMode(int mode)
{
    this->sampleDistribution.setDistributionMode(mode);
}

float GenotypeResult::getQuality()
{
    return this->quality;
}

void GenotypeResult::calculateReadStats()
{
    std::unordered_map<std::string, float> groupQualities;
    for (auto & key : this->readPairGroups) 
    {
        std::vector<float> qualities;
        for (int & i : key.second)
            qualities.push_back(this->mappingQualities[i]);
        groupQualities[key.first] = mean(qualities);
    }
    this->groupQualities = groupQualities;
}

std::unordered_map<std::string, float> GenotypeResult::getReadStats()
{
    return this->groupQualities;
}

void GenotypeResult::clearData()
{
    this->sampleDistribution = GenotypeDistribution();
    std::vector<int>().swap(this->mappingQualities);
    std::unordered_map<std::string, std::vector<int>>().swap(this->readPairGroups);
    std::vector<float>().swap(this->templateWeights);
    std::vector<float>().swap(this->bootstrappedQualities);
    std::vector<std::vector<float>>().swap(this->bootstrappedLikelihoods);
    std::vector<std::vector<float>>().swap(this->templateProbabilities);
    std::vector<std::string>().swap(this->genotypeNames);
    std::vector<float>().swap(genotypeLikelihoods);
    std::vector<float>().swap(genotypeLikelihoods);
}