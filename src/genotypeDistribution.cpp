#include "genotypeDistribution.hpp"
#include "insertSizeDistribution.hpp"
#include <cstdint>

GenotypeDistribution::GenotypeDistribution()
{
    this->normalizationFactor = 1;
    this->minProbability = 0.0;
}

GenotypeDistribution::GenotypeDistribution(const Eigen::SparseMatrix<float, Eigen::RowMajor, int64_t> & dist, std::unordered_map<std::string, int> & groupIndices, int64_t sMin, int64_t sMax)
{
    this->normalizationFactor = 1;
    this->minProbability = 0.0;

    for (auto & g : groupIndices)
    {
        if (dist.row(g.second).nonZeros() == 0)
            continue;
        this->distributions[g.first] = InsertSizeDistribution(sMin, sMax, dist.row(g.second));
        this->distributions[g.first].smoothDistribution(5);
    }
        
    scaleDistribution();
    for (auto & dist : this->distributions)
        this->distributionProbabilities[dist.first] = dist.second.getIntegral();
    
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

    scaleDistribution();
    return *this;
}

GenotypeDistribution & GenotypeDistribution::operator*=(float & factor)
{
    for (auto it : this->distributions)
        this->distributions[it.first] *= factor;
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


void GenotypeDistribution::addInsertSizeProbability(int64_t insertSize, std::string orientation, std::string junctionString, std::string breakpointString, std::string bridgeString, float probability)
{
    std::string group = GenotypeDistribution::determineGroup(
        orientation, junctionString, breakpointString, bridgeString
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


void GenotypeDistribution::calculateNormalizationFactor()
{
    float dSum = 0;
    for (auto it : this->distributions) 
    {
        float temp = it.second.getIntegral();
        this->distributionProbabilities[it.first] = temp;
        dSum += temp;
    }
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

float GenotypeDistribution::getProbability(int64_t insertSize, std::string orientation, std::string junctionString, std::string breakpointString, std::string bridgeString, bool & outlier)
{

    std::string group = GenotypeDistribution::determineGroup(
        orientation, junctionString, breakpointString, bridgeString
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
    int64_t bounds = 20000;
    f.open(prefix.append(".txt"));
    if (f.is_open())
    {
        for (auto & it : this->distributions) {
            float oob_low = 0;
            float oob_high = 0;

            for (Eigen::SparseVector<float, 0, int64_t>::InnerIterator it1(it.second.getDistributionVector(), 0); it1; ++it1)
            {
                float p = it1.value();
                int64_t idx = it1.row();
                if (idx + it.second.getMinInsertSize() < -bounds)
                    oob_low += p;
                else if (idx + it.second.getMinInsertSize() > bounds)
                    oob_high += p;
                else
                    f << idx + it.second.getMinInsertSize() << "\t" << p << "\t" << it.first << std::endl; 
            }

            if (oob_low > 0)
                f << 0 << "\t" << oob_low << "\t" << it.first << "_oob_low" << std::endl;
            if (oob_high > 0)
                f << 0 << "\t" << oob_high << "\t" << it.first << "_oob_high" << std::endl;
        }
    }
    f.close();
}

void GenotypeDistribution::writeDistributionBinned(std::string prefix)
{
    int binSize = 20;
    int bounds = 20000;

    std::ofstream f;
    f.open(prefix.append(".txt"));
    if (f.is_open())
    {
        // for each distribution: 
        for (auto & it : this->distributions) {
            int64_t maxIS = it.second.getMaxInsertSize();
            int64_t minIS = it.second.getMinInsertSize();

            float oob_low = 0;
            float oob_high = 0;

            //  create a binned vector of probabilities
            Eigen::SparseVector<float, 0, int64_t> binnedDist;

            int64_t npBins = maxIS / binSize;
            if (maxIS % binSize)
                npBins++;
            int64_t nnBins = 0;
            if (minIS < 0) {
                nnBins = (-minIS) / binSize;
                if (minIS % binSize)
                    nnBins++;
            }
            binnedDist.resize(nnBins + npBins);
            binnedDist.reserve(it.second.getDistributionVector().nonZeros() / binSize + 1);

            for (Eigen::SparseVector<float, 0, int64_t>::InnerIterator it1(it.second.getDistributionVector(), 0); it1; ++it1)
            {
                float p = it1.value();
                int64_t idx = it1.row();
                int64_t newIdx = idx / binSize;
                binnedDist.coeffRef(newIdx) += p;
            }

            //  write as above
            for (Eigen::SparseVector<float, 0, int64_t>::InnerIterator it1(binnedDist, 0); it1; ++it1)
            { 
                float p = it1.value();
                int64_t idx = it1.row();
                int64_t sBin = (idx - nnBins) * binSize;

                if (sBin < -bounds)
                    oob_low += p;
                else if (sBin > bounds)
                    oob_high += p;
                else
                    f << sBin << "\t" << p << "\t" << it.first << std::endl; 
            }

            if (oob_low > 0)
                f << 0 << "\t" << oob_low << "\t" << it.first << "_oob_low" << std::endl;
            if (oob_high > 0)
                f << 0 << "\t" << oob_high << "\t" << it.first << "_oob_high" << std::endl;
        }
    }
    f.close();
}

void GenotypeDistribution::scaleDistribution()
{
    calculateNormalizationFactor();
    for (auto it : this->distributions)
        this->distributions[it.first] *= this->normalizationFactor;
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

int64_t GenotypeDistribution::getMinInsertSize()
{
	int64_t minIS = 0;
	for (auto it : this->distributions)
		minIS = std::min(minIS, it.second.getMinInsertSize());
	return minIS;
}

int64_t GenotypeDistribution::getMaxInsertSize()
{
	int64_t maxIS = 0;
	for (auto it : this->distributions)
		maxIS = std::max(maxIS, it.second.getMaxInsertSize());
	return maxIS;
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

    for (auto it : this->distributions)
    {
        for (Eigen::SparseVector<float, 0, int64_t>::InnerIterator it1(it.second.getDistributionVector(), 0); it1; ++it1)
        {
            if (it1.value() <= this->minProbability)
                continue;
            
            float x1 = it1.value();
            float x2 = std::log10(x1);
            float x3 = std::log10(std::max(this->minProbability, other.getDistributions()[it.first].getInsertSizeProbability(it1.row() + it.second.getMinInsertSize())));
            kld += x1 * (x2 - x3);
        }
    }

    kld *= 10; //phred scaling
    return kld;
}
