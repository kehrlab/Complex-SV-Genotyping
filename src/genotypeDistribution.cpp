#include "genotypeDistribution.hpp"
#include "insertSizeDistribution.hpp"

GenotypeDistribution::GenotypeDistribution()
{
    this->normalizationFactor = 1;
    this->minProbability = 0.0;
}

GenotypeDistribution::GenotypeDistribution(const Eigen::SparseMatrix<float> & dist, std::unordered_map<std::string, int> & groupIndices, int sMin, int sMax)
{
    this->normalizationFactor = 1;
    this->minProbability = 0.0;

    for (auto & g : groupIndices)
    {
        int nonZero = 0;
        std::vector<float> temp(dist.rows(), 0);
        for (int i = 0; i < dist.rows(); ++i)
        {
            float p = dist.coeff(i, g.second);
            if (p > 0)
            {
                ++nonZero;
                temp[i] = p;
            }
        }

        if (nonZero == 1 && temp[0 - sMin] > 0)
            this->interChromosomeCounts[g.first] = temp[0 - sMin];
        else
            this->distributions[g.first] = InsertSizeDistribution(sMin, sMax, temp);
    }
    
    for (auto & dist : this->distributions)
        this->distributionProbabilities[dist.first] = dist.second.getIntegral();
    this->distributionProbabilities["interchromosome"] = getTotalInterChromosomeProbability();
    
    calculateMinProbability();
}

GenotypeDistribution & GenotypeDistribution::operator+=(GenotypeDistribution & rhs)
{	
    for (auto it : rhs.getDistributions()) {
        if (this->distributions.find(it.first) == this->distributions.end())
            this->distributions[it.first] = it.second;
        else
            this->distributions[it.first] += it.second;
    }

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

void GenotypeDistribution::addInsertSizeProbability(int insertSize, std::string orientation, std::string junctionString, std::string breakpointString, std::string chromosomeString, float probability)
{
    if (chromosomeString != "")
    {
        addInterChromosomeProbability(chromosomeString, probability);
        return;
    }

    std::string group = GenotypeDistribution::determineGroup(
        orientation, junctionString, breakpointString, chromosomeString
        );

    if (group == "")
    {
        std::cerr << "Could not assign group to read pair" << std::endl;
        return;
    }

    if (this->distributions.find(group) == this->distributions.end())
        this->distributions[group] = InsertSizeDistribution();
    this->distributions[group].addInsertSizeProbability(insertSize, probability);
    return;
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

float GenotypeDistribution::getProbability(int insertSize, std::string orientation, std::string junctionString, std::string breakpointString, std::string chromosomeString, bool & outlier)
{
    if (chromosomeString != "")
        return getInterChromosomeProbability(chromosomeString);

    std::string group = GenotypeDistribution::determineGroup(
        orientation, junctionString, breakpointString, chromosomeString
        );

    if (this->distributions.find(group) != this->distributions.end())
        return std::max(this->distributions[group].getInsertSizeProbability(insertSize), this->minProbability);

    // otherwise
    outlier = true;
    return this->minProbability;
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
    {
        return this->interChromosomeCounts[identifier];
    } else
    {
        std::cerr << "Did not find inter-chromosome group " << identifier << std::endl;
        return this->minProbability;
    }
}

std::unordered_map<std::string, float> GenotypeDistribution::getInterChromosomeProbabilites()
{
    return this->interChromosomeCounts;
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
