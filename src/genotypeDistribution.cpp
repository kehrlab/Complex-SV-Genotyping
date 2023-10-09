#include "genotypeDistribution.hpp"
#include "insertSizeDistribution.hpp"
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

GenotypeDistribution::GenotypeDistribution()
{
    this->normalizationFactor = 1;
    this->minProbability = 0.0;
    this->mode = 2;
}

GenotypeDistribution::GenotypeDistribution(std::vector<std::string> contigNames, int mode)
{
    this->normalizationFactor = 1;
    this->minProbability = 0.0;
    this->mode = 2;
    setPossibleContigs(contigNames);
    setDistributionMode(mode);
}

GenotypeDistribution::GenotypeDistribution(const Eigen::SparseMatrix<float> & dist, std::unordered_map<std::string, int> & groupIndices, int sMin, int sMax)
{
    this->normalizationFactor = 1;
    this->minProbability = 0.0;
    this->mode = 2;

    for (auto & g : groupIndices)
    {
        std::vector<float> temp(dist.rows(), 0);
        for (int i = 0; i < dist.rows(); ++i)
            temp[i] = dist.coeff(i, g.second);
        this->distributions[g.first] = InsertSizeDistribution(sMin, sMax, temp);
    }
    calculateMinProbability();
}

void GenotypeDistribution::setDistributionMode(int mode)
{
    if (mode >= 0 && mode <=3)
        this->mode = mode;
    else {
        std::cout << "Invalid mode. Set to default (2)." << std::endl;
    }
    this->distributions.erase(this->distributions.begin(), this->distributions.end());
    this->distributionProbabilities.erase(this->distributionProbabilities.begin(), this->distributionProbabilities.end());
}

GenotypeDistribution & GenotypeDistribution::operator+=(GenotypeDistribution & rhs)
{	
    for (auto it : rhs.getDistributions()) {
        if (this->distributions.find(it.first) == this->distributions.end())
            this->distributions[it.first] = it.second;
        else
            this->distributions[it.first] += it.second;
    }

    std::unordered_set<std::string> rhs_contigs = rhs.getPossibleContigs();
    for (auto cName : rhs_contigs)
        this->contigs.insert(cName);

    std::unordered_map<std::string, float> rhs_probs = rhs.getInterChromosomeProbabilites();
    for (auto it : rhs_probs)
    {
        if (this->interChromosomeCounts.find(it.first) == this->interChromosomeCounts.end())
            this->interChromosomeCounts[it.first] = 0.0;
        this->interChromosomeCounts[it.first] += it.second;
    }

    scaleDistribution();
    return *this;
}

GenotypeDistribution & GenotypeDistribution::operator*=(float & factor)
{
    for (auto it : this->distributions)
        this->distributions[it.first] *= factor;
    for (auto it : this->interChromosomeCounts)
        this->interChromosomeCounts[it.first] *= factor;
    return *this;
}

GenotypeDistribution operator+(GenotypeDistribution lhs, GenotypeDistribution rhs)
{
    lhs += rhs;
    return lhs;
}

GenotypeDistribution operator*(GenotypeDistribution lhs, float factor)
{
    lhs *= factor;
    return lhs;
}

GenotypeDistribution operator*(float factor, GenotypeDistribution rhs)
{
    rhs *= factor;
    return rhs;
}

void GenotypeDistribution::addInsertSizeProbability(int insertSize, std::string orientation, bool split, bool spanning, bool interChromosome, std::string regionString, std::string junctionString, std::string breakpointString, std::vector<std::vector<std::string>> rNamePairs, float probability)
{
    if (interChromosome)
    {
        if (rNamePairs.size() < 1)
            std::cout << "WARNING: No chromosome names for inter-chromosome read pair" << std::endl;
        for (unsigned i = 0; i < rNamePairs.size(); ++i)
        {
            if (rNamePairs[i].size() != 2)
                throw std::runtime_error("Invalid number of chromosome names for inter-chromosome read pair");
            std::string contigIdentifier = getContigIdentifier(rNamePairs[i][0], rNamePairs[i][1]);
            if (contigIdentifier != "")
                addInterChromosomeProbability(contigIdentifier, probability);
        }
        if (insertSize == 0)
            return;
    }

    std::vector<std::string> distKeys = GenotypeDistribution::determineDistributionKeys(orientation, split, spanning, this->mode, regionString, junctionString, breakpointString);
    for (std::string distKey : distKeys) {
        if (this->distributions.find(distKey) == this->distributions.end()) {
            this->distributions[distKey] = InsertSizeDistribution();
            this->distributionProbabilities[distKey] = 0.;
        }
        this->distributions[distKey].addInsertSizeProbability(insertSize, probability);
    }
    return;
}

std::vector<std::string> GenotypeDistribution::determineDistributionKeys(std::string orientation, bool split, bool spanning, int mode, std::string regionString, std::string junctionString, std::string breakpointString)
{
    std::vector<std::string> distKeys;

    if (mode == 1) {
        if (orientation == "FR" || orientation == "RF")
        {
            if (split)
                distKeys.push_back("splitRF");
            if (spanning)
                distKeys.push_back("spanningRF");
            if (!split && !spanning)
                distKeys.push_back("RF");
        }
        else if (orientation == "RR")
        {
            if (split)
                distKeys.push_back("splitRR");
            if (spanning)
                distKeys.push_back("spanningRR");
            if (!split && !spanning)
                distKeys.push_back("RR");
        } 
        else if (orientation == "FF")
        {
            if (split)
                distKeys.push_back("splitFF");
            if (spanning)
                distKeys.push_back("spanningFF");
            if (!split && !spanning)
                distKeys.push_back("FF");
        }
    } else {
        if (split)
            distKeys.push_back("split");
        else if (spanning)
            distKeys.push_back("spanning");
        else {
            if (orientation == "FR" || orientation == "RF")
                distKeys.push_back("RF");
            else
                distKeys.push_back(orientation);
        }
    }

    if (mode == 2)
    {
        for (int i = 0; i < distKeys.size(); ++i)
        {
            if (distKeys[i] == "spanning" && breakpointString != "")
                distKeys[i] = distKeys[i] + "_" + breakpointString;
            else if (distKeys[i] == "split" && junctionString != "")
                distKeys[i] = distKeys[i] + "_" + junctionString;
        }
    }
    
    return distKeys;
}

float GenotypeDistribution::getTotalInterChromosomeProbability()
{
    float total = 0.0;
    for (auto it : this->interChromosomeCounts)
        total += it.second;
    return total;
}

void GenotypeDistribution::calculateNormalizationFactor()
{
    float dSum = 0;
    for (auto it : this->distributions) 
    {
        float temp = it.second.getIntegral();
        this->distributionProbabilities[it.first] = temp;
        dSum += temp;
    }
    dSum += getTotalInterChromosomeProbability();
    this->normalizationFactor = 1 / dSum;

    for (auto it : this->distributionProbabilities)
        this->distributionProbabilities[it.first] = this->distributionProbabilities[it.first] * this->normalizationFactor;
    return;
}

float GenotypeDistribution::getNormalizationFactor()
{
    calculateNormalizationFactor();
    return this->normalizationFactor;
}

float GenotypeDistribution::getProbability(int insertSize, std::string orientation, bool split, bool spanning, bool interChromosome, std::string regionString, std::string junctionString, std::string breakpointString, std::vector<std::vector<std::string>> rNamePairs, bool useInsertSize)
{
    if (interChromosome)
    {
        float prob = this->minProbability;
        for (unsigned i = 0; i < rNamePairs.size(); ++i)
        {
            if (rNamePairs[i].size() != 2)
                throw std::runtime_error("Invalid number of chromosome names for inter-chromosome read pair.");
            prob = std::max(prob, getInterChromosomeProbability(getContigIdentifier(rNamePairs[i][0], rNamePairs[i][1])));
        }
        return prob;
    }

    std::vector<std::string> distKeys = determineDistributionKeys(orientation, split, spanning, this->mode, regionString, junctionString, breakpointString);

    float p = 0.0;
    if (useInsertSize)
    {
        if (this->mode != 3)
        {
            if (distKeys.size() == 2)
                p = std::max(this->distributions[distKeys[0]].getInsertSizeProbability(insertSize), this->distributions[distKeys[1]].getInsertSizeProbability(insertSize));
            else if (distKeys.size() == 1)
                p = this->distributions[distKeys[0]].getInsertSizeProbability(insertSize);
        } else {
            if (distKeys.size() == 2)
            {
                p = std::max(this->distributions[distKeys[0]].getInsertSizeProbability(insertSize), this->distributions[distKeys[1]].getInsertSizeProbability(insertSize));
            }
            else if (distKeys.size() == 1)
            {
                if (distKeys[0] == "RF")
                    p = this->distributions[distKeys[0]].getInsertSizeProbability(insertSize);
                else
                    p = this->distributionProbabilities[distKeys[0]];
            }
        }
    } else {
        if (distKeys.size() == 2)
            p = std::max(this->distributionProbabilities[distKeys[0]], this->distributionProbabilities[distKeys[1]]);
        else if (distKeys.size() == 1)
            p = this->distributionProbabilities[distKeys[0]];
    }
    
    return std::max(this->minProbability, p);
}

std::unordered_map<std::string, InsertSizeDistribution> & GenotypeDistribution::getDistributions()
{
    return this->distributions;
}

void GenotypeDistribution::writeDistribution(std::string prefix)
{
    std::ofstream f;
    f.open(prefix.append(".txt"));
    if (f.is_open())
    {
        for (auto it : this->distributions)
            for (int i = it.second.getMinInsertSize(); i <= it.second.getMaxInsertSize(); ++i)
                f << i << "\t" << it.second.getInsertSizeProbability(i) << "\t" << it.first << std::endl;
        for (auto it : this->interChromosomeCounts)
            f << 0 << "\t" << it.second << "\t" << it.first << std::endl;
    }
    f.close();
}

void GenotypeDistribution::writeDistributionBinned(std::string prefix)
{
    int binSize = 20;
    int minIS = this->distributions["RF"].getMinInsertSize();
    int maxIS = 0;
    for (auto it : this->distributions)
        maxIS = std::max(maxIS, it.second.getMaxInsertSize());
    std::ofstream f;
    f.open(prefix.append(".txt"));
    if (f.is_open())
    {
        std::unordered_map<std::string, float> distProbs;
        for (auto it : this->distributions)
            distProbs[it.first] = 0.0;

        int j = 0;
        int iS = minIS;
        for (int i = minIS; i < maxIS; ++i)
        {
		
            ++j;
            if (j % binSize == 0)
            {
                for (auto it : distProbs) {
                    f << i << "\t" << it.second << "\t" << it.first << std::endl;
                    distProbs[it.first] = 0.0;
                }
                iS = i;
            }
            for (auto it : this->distributions)
                distProbs[it.first] += it.second.getInsertSizeProbability(i);
        }
        for (auto it : this->interChromosomeCounts)
            f << 0 << "\t" << it.second << "\t" << it.first << std::endl;
    }
    f.close();
}

void GenotypeDistribution::scaleDistribution()
{
    calculateNormalizationFactor();
    for (auto it : this->distributions)
        this->distributions[it.first] *= this->normalizationFactor;
    for (auto it : this->interChromosomeCounts)
	    this->interChromosomeCounts[it.first] = this->interChromosomeCounts[it.first] * this->normalizationFactor;
    this->normalizationFactor = 1.0;
    calculateMinProbability();
}

void GenotypeDistribution::smoothDistribution()
{
    for (auto it : this->distributions) {
        this->distributions[it.first].smoothDistribution(5);
        this->distributions[it.first].smoothDistribution(10);
        this->distributions[it.first].smoothDistribution(5);
    }
}

int GenotypeDistribution::getMinInsertSize()
{
	int minIS = 0;
	for (auto it : this->distributions)
		minIS = std::min(minIS, it.second.getMinInsertSize());
	return minIS;
}

int GenotypeDistribution::getMaxInsertSize()
{
	int maxIS = 0;
	for (auto it : this->distributions)
		maxIS = std::max(maxIS, it.second.getMaxInsertSize());
	return maxIS;
}

void GenotypeDistribution::setPossibleContigs(std::vector<std::string> cNames)
{
    this->contigs.erase(this->contigs.begin(), this->contigs.end());
    for (auto cName : cNames)
        this->contigs.insert(cName);
}

std::string GenotypeDistribution::getContigIdentifier(std::string rName1, std::string rName2)
{
    std::string identifier = "";

    if (this->contigs.size() > 0)
        if (this->contigs.find(rName1) == this->contigs.end() || this->contigs.find(rName2) == this->contigs.end())
            return "";

    if (rName1 <= rName2)
        identifier = rName1 + rName2;
    else
        identifier = rName2 + rName1;
    return identifier;
}

void GenotypeDistribution::addInterChromosomeProbability(std::string identifier, float probability)
{
    auto it = this->interChromosomeCounts.find(identifier);
    if (it == interChromosomeCounts.end())
        this->interChromosomeCounts[identifier] = 0.0;
    this->interChromosomeCounts[identifier] += probability;
}

float GenotypeDistribution::getInterChromosomeProbability(std::string identifier)
{
    auto it = this->interChromosomeCounts.find(identifier);
    if (it != this->interChromosomeCounts.end())
        return this->interChromosomeCounts[identifier];
    else
        return this->minProbability;
}

std::unordered_map<std::string, float> GenotypeDistribution::getInterChromosomeProbabilites()
{
    return this->interChromosomeCounts;
}

std::unordered_set<std::string> GenotypeDistribution::getPossibleContigs()
{
    return this->contigs;
}

void GenotypeDistribution::calculateMinProbability()
{
    for (auto it : this->distributions)
        this->distributions[it.first].findMinProbability();

    this->minProbability = 0.001;
    for (auto it : this->distributions)
        if (it.second.getMinProbability() > 0)
            this->minProbability = std::min(it.second.getMinProbability(), this->minProbability);
    this->minProbability /= 2;
}

float GenotypeDistribution::getMinProbability()
{
    return this->minProbability;
}

void GenotypeDistribution::setMinProbability(float minProbability)
{
    this->minProbability = minProbability;
}

float GenotypeDistribution::calculateKLD(GenotypeDistribution & other)
{
    float kld = 0.0;
    // problem: values with p = 0 do not count -> REF/REF likely
    for (auto it : this->distributions)
    {
        int min = it.second.getMinInsertSize();
        int max = it.second.getMaxInsertSize();
        for (int x = min; x <= max; ++x) {
            if (it.second.getInsertSizeProbability(x) <= this->minProbability)
                continue;
            float x1 = it.second.getInsertSizeProbability(x);
            float x2 = std::log10(x1);
            float x3 = std::log10(std::max(other.getDistributions()[it.first].getInsertSizeProbability(x), this->minProbability));
            kld += x1 * (x2 - x3);
        }
    }
    for (auto it : this->interChromosomeCounts)
    {
        float x1 = it.second;
        float x2 = std::log10(x1);
        float x3 = std::log10(other.getInterChromosomeProbability(it.first));
        kld += x1 * (x2-x3);
    }
    kld *= 10; //phred scaling
    return kld;
}
