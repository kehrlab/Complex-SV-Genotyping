#include "insertSizeDistribution.hpp"
#include <fstream>

InsertSizeDistribution::InsertSizeDistribution()
{
    this->maxInsertSize = 1000;
    this->minInsertSize = 1;
    this->minValue = 0.0;
    std::vector<float> distribution(1000, 0.0);
    this->distribution = distribution;
}

InsertSizeDistribution::InsertSizeDistribution(std::vector<int> insertSizes)
{
    this->minInsertSize = 1000;
    this->maxInsertSize = 1;
    this->minValue = 0.0;
    if (insertSizes.size() > 0) 
    {
        inferMinAndMaxValues(insertSizes);
        std::vector<float> distribution = std::vector<float>(this->maxInsertSize - this->minInsertSize + 1, 0.0);
        this->distribution = distribution;
        calculateDistributionMean(insertSizes);
        calculateDistributionSD(insertSizes);
        generateDistributionFromInsertSizes(insertSizes);
    }
    else 
    {
        std::vector<float> distribution(1000, 0.0);
        this->distribution = distribution;
    }
}

InsertSizeDistribution::InsertSizeDistribution(int min, int max, std::vector<float> distribution)
{
    this->minValue = 0.0;
    this->minInsertSize = min;
    this->maxInsertSize = max;
    this->distribution = distribution;
}

InsertSizeDistribution & InsertSizeDistribution::operator+=(InsertSizeDistribution & rhs)
{
    for (int i = rhs.getMinInsertSize(); i <= rhs.getMaxInsertSize(); ++i)
        addInsertSizeProbability(i, rhs.getInsertSizeProbability(i));
    return *this;
}

InsertSizeDistribution operator+(InsertSizeDistribution lhs, InsertSizeDistribution rhs)
{
    lhs += rhs;
    return lhs;
}

InsertSizeDistribution & InsertSizeDistribution::operator*=(float & factor)
{
    for (unsigned i = 0; i < this->distribution.size(); ++i)
        this->distribution[i] *= factor;
    return *this;
}

InsertSizeDistribution operator*(float factor, InsertSizeDistribution dist)
{
    dist *= factor;
    return dist;
}

InsertSizeDistribution operator*(InsertSizeDistribution dist, float factor)
{
    dist *= factor;
    return dist;
}

void InsertSizeDistribution::addInsertSizeProbability(int insertSize, float p)
{
    if (insertSize >= this->minInsertSize && insertSize <= this->maxInsertSize)
    {
        this->distribution[insertSize - this->minInsertSize] += p;
    } else if (this->minInsertSize > insertSize)
    {
        std::vector<float> tempDistribution(this->maxInsertSize - insertSize + 101, 0);
        for (unsigned i = 0; i < this->distribution.size(); ++i)
        {
            tempDistribution[this->minInsertSize-(insertSize-100) +i] = this->distribution[i];
        }
        tempDistribution[100] = p;

        this->minInsertSize = insertSize - 100;
        this->distribution = tempDistribution;
    } else if (this->maxInsertSize < insertSize)
    {
        std::vector<float> tempDistribution(insertSize + 101 - this->minInsertSize, 0);
        for (unsigned i = 0; i < this->distribution.size(); ++i)
        {
            tempDistribution[i] = this->distribution[i];
        }
        tempDistribution[insertSize - this->minInsertSize] = p;

        this->maxInsertSize = insertSize + 100;
        this->distribution = tempDistribution;
    }
}

void InsertSizeDistribution::inferMinAndMaxValues(std::vector<int> const & insertSizes)
{
    for (auto it : insertSizes)
    {
        if (it < this->minInsertSize)
            this->minInsertSize = it;
        if (it > this->maxInsertSize)
            this->maxInsertSize = it;
    }
    this->minInsertSize = std::max(this->minInsertSize - 50, 1);
    this->maxInsertSize += 50;
}

void InsertSizeDistribution::generateDistributionFromInsertSizes(std::vector<int> const & insertSizes)
{
    distribution = std::vector<float>(this->maxInsertSize - this->minInsertSize + 1, 0);
    for (int insertSize = this->minInsertSize; insertSize <= this->maxInsertSize; ++insertSize)
    {
        this->distribution[insertSize - this->minInsertSize] = getNumberOfInsertSizes(insertSizes, insertSize);
        if (this->distribution[insertSize - this->minInsertSize] == 0)
            this->distribution[insertSize - this->minInsertSize] = this->minValue;
    }
    smoothDistribution(21);
    smoothDistribution(41);
    smoothDistribution(11);
    scaleDistribution();
}

int InsertSizeDistribution::getNumberOfInsertSizes(std::vector<int> const & insertSizes, int insertSize)
{
    int occurrences = 0;
    for (auto it : insertSizes)
        if (it == insertSize)
            ++occurrences;
    return occurrences;
}

void InsertSizeDistribution::calculateDistributionMean(std::vector<int> const & insertSizes)
{
    this->insertMean = 0;
    for (auto it : insertSizes)
        this->insertMean += it;
    this->insertMean /= insertSizes.size();
}

void InsertSizeDistribution::calculateDistributionSD(std::vector<int> const & insertSizes)
{
    this->insertSD = 0;
    for (auto it : insertSizes)
        this->insertSD += (it - this->insertMean) * (it-this->insertMean);
    this->insertSD = std::sqrt(this->insertSD / insertSizes.size());
}

void InsertSizeDistribution::scaleDistribution()
{
    float distributionSum = 0;
    for (auto it : this->distribution)
        distributionSum += it;
    for (int i = 0; i < this->distribution.size(); ++i)
        this->distribution[i] /= distributionSum;
    findMinProbability();
}

void InsertSizeDistribution::smoothDistribution(int kernelSize)
{
    float *kernel = new float[kernelSize];
    std::vector<float> smoothed_values(this->distribution.size(), 0);

    for (int i = - (kernelSize-1)/2; i <= (kernelSize-1)/2; ++i)
    {
        kernel[i+(kernelSize-1)/2] = 1 / std::sqrt(2 * M_PI) * std::exp(- (i*i)/((kernelSize-1)));
    }

    for (unsigned i = 0; i < this->distribution.size(); ++i)
    {
        for (int j = 0; j < kernelSize; ++j)
        {
            if (i+j-(kernelSize-1)/2 < 0 || i+j-(kernelSize-1)/2 >= this->distribution.size())
                continue; 
            smoothed_values[i] += kernel[j] * this->distribution[i+j-(kernelSize-1)/2];
        }
        if (smoothed_values[i] == 0)
            smoothed_values[i] = this->minValue;
    }
    this->distribution = smoothed_values;

    delete[] kernel;
}

void InsertSizeDistribution::writeDistribution(std::string filename)
{
    if (filename == "") {
        for (int i = this->minInsertSize; i <= maxInsertSize; ++i)
            std::cout << i << "\t" << this->distribution[i-this->minInsertSize] << std::endl;
    } else {
        std::ofstream f;
        f.open(filename);
        if (f.is_open())
        {
            for (int i = this->minInsertSize; i <= maxInsertSize; ++i)
                f << i << "\t" << this->distribution[i-this->minInsertSize] << std::endl;
        }
        f.close();
    }
}

double InsertSizeDistribution::getInsertMean() const
{
    return this->insertMean;
}

double InsertSizeDistribution::getInsertSD() const
{
    return this->insertSD;
}

float InsertSizeDistribution::getInsertSizeProbability(int insertSize)
{
    if (insertSize >= this->minInsertSize && insertSize <= this->maxInsertSize)
        return this->distribution[insertSize - this->minInsertSize];
    else
        return 0;
}

int InsertSizeDistribution::getMinInsertSize()
{
    return this->minInsertSize;
}

int InsertSizeDistribution::getMaxInsertSize()
{
    return this->maxInsertSize;
}

float InsertSizeDistribution::getIntegral()
{
    float sum = 0;
    for (auto it : this->distribution)
        sum += it;
    return sum;
}

void InsertSizeDistribution::findMinProbability()
{
    this->minValue = 1.0;
    for (float p : this->distribution)
        if (p > 0 && p < this->minValue)
            this->minValue = p;
    if (this->minValue == 1.0)
        this->minValue = 0.0;
}

float InsertSizeDistribution::getMinProbability()
{
    return this->minValue;
}