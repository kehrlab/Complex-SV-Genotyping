#include "insertSizeDistribution.hpp"

InsertSizeDistribution::InsertSizeDistribution()
{
    this->maxInsertSize = 1000;
    this->minInsertSize = 1;
    this->minValue = 0.0;

    this->distribution = Eigen::SparseVector<float, 0, int64_t>(1000);
}

InsertSizeDistribution::InsertSizeDistribution(std::vector<int64_t> insertSizes)
{
    this->minInsertSize = 1000;
    this->maxInsertSize = 1;
    this->minValue = 0.0;
    if (insertSizes.size() > 0) 
    {
        inferMinAndMaxValues(insertSizes);
        this->distribution = Eigen::SparseVector<float, 0, int64_t>(this->maxInsertSize - this->minInsertSize + 1);
        generateDistributionFromInsertSizes(insertSizes); // this also calculates mean insert size
        calculateDistributionSD(insertSizes);
    }
    else 
    {
        this->distribution = Eigen::SparseVector<float, 0, int64_t>(1000);
    }
}

InsertSizeDistribution::InsertSizeDistribution(int64_t min, int64_t max, Eigen::SparseVector<float, 0, int64_t> distribution)
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
    this->distribution *= factor;
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

void InsertSizeDistribution::addInsertSizeProbability(int64_t insertSize, float p)
{
    // resize distribution vector if required
    if (insertSize < this->minInsertSize || insertSize > this->maxInsertSize)
    {
        // determine new borders and resize
        int64_t newMin = std::min(this->minInsertSize, insertSize - 100);
        int64_t newMax = std::max(this->maxInsertSize, insertSize + 100);
        this->distribution.conservativeResize(newMax - newMin + 1);

        // shift all entries if resize was at the beginning
        if (newMin < this->minInsertSize)
        {
            int64_t diff = this->minInsertSize - newMin;
            for (int64_t i = 0; i < this->distribution.nonZeros(); ++i)
                *(this->distribution.innerIndexPtr() + i) += diff;
            this->minInsertSize = newMin;
        } else {
            this->maxInsertSize = newMax;
        }
    } 
    // add probability
    distribution.coeffRef(insertSize - this->minInsertSize) += p;
}

void InsertSizeDistribution::inferMinAndMaxValues(std::vector<int64_t> const & insertSizes)
{
    for (auto it : insertSizes)
    {
        if (it < this->minInsertSize)
            this->minInsertSize = it;
        if (it > this->maxInsertSize)
            this->maxInsertSize = it;
    }
    this->minInsertSize = std::max(this->minInsertSize - 100, (int64_t) 1);
    this->maxInsertSize += 100;
}

void InsertSizeDistribution::generateDistributionFromInsertSizes(std::vector<int64_t> const & insertSizes)
{
    distribution = Eigen::SparseVector<float, 0, int64_t>(this->maxInsertSize - this->minInsertSize + 1);
    this->insertMean = 0;    
    for (uint32_t i = 0; i < insertSizes.size(); ++i) {
        this->distribution.coeffRef(insertSizes[i] - this->minInsertSize) += 1;
        this->insertMean += insertSizes[i];
    }
    this->insertMean /= insertSizes.size();
    
    smoothDistribution(21);
    smoothDistribution(41);
    smoothDistribution(11);
    scaleDistribution();
}

void InsertSizeDistribution::calculateDistributionMean(std::vector<int64_t> const & insertSizes)
{
    this->insertMean = 0;
    for (auto it : insertSizes)
        this->insertMean += it;
    this->insertMean /= insertSizes.size();
}

void InsertSizeDistribution::calculateDistributionSD(std::vector<int64_t> const & insertSizes)
{
    this->insertSD = 0;
    for (auto it : insertSizes)
        this->insertSD += (it - this->insertMean) * (it-this->insertMean);
    this->insertSD = std::sqrt(this->insertSD / insertSizes.size());
}

void InsertSizeDistribution::scaleDistribution()
{
    float distributionSum = this->distribution.sum();
    this->distribution /= distributionSum;
    findMinProbability();
}

void InsertSizeDistribution::smoothDistribution(int kernelSize)
{
    Eigen::SparseVector<float, 0, int64_t> smoothed_values(this->distribution.rows());
    smoothed_values.reserve(this->distribution.nonZeros() * 5);

    std::vector<float> kernel(5, 0);
    for (int i = - (kernelSize-1)/2; i <= (kernelSize-1)/2; ++i)
        kernel[i+(kernelSize-1)/2] = 1 / std::sqrt(2 * M_PI) * std::exp(- (i*i)/((kernelSize-1)));

    for (Eigen::SparseVector<float, 0, int64_t>::InnerIterator it(this->distribution, 0); it; ++it)
    {
        int64_t idx = it.row();
        float value = 0;
        for (int j = 0; j < kernelSize; ++j)
        {
            if (idx+j-(kernelSize-1)/2 < 0 || idx+j-(kernelSize-1)/2 >= this->distribution.size())
                continue; 
            value += kernel[j] * this->distribution.coeff(idx+j-(kernelSize-1)/2);
        }
        smoothed_values.insert(idx) = value;
    }
    this->distribution = smoothed_values;
}

void InsertSizeDistribution::writeDistribution(std::string filename)
{
    if (filename == "") {
        for (int i = this->minInsertSize; i <= maxInsertSize; ++i)
            std::cout << i << "\t" << this->distribution.coeff(i-this->minInsertSize) << std::endl;
    } else {
        std::ofstream f;
        f.open(filename);
        if (f.is_open())
        {
            for (int i = this->minInsertSize; i <= maxInsertSize; ++i)
                f << i << "\t" << this->distribution.coeff(i - this->minInsertSize) << std::endl;
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

float InsertSizeDistribution::getInsertSizeProbability(int64_t insertSize)
{
    if (insertSize >= this->minInsertSize && insertSize <= this->maxInsertSize)
        return this->distribution.coeff(insertSize - this->minInsertSize);
    else
        return 0;
}

int64_t InsertSizeDistribution::getMinInsertSize()
{
    return this->minInsertSize;
}

int64_t InsertSizeDistribution::getMaxInsertSize()
{
    return this->maxInsertSize;
}

float InsertSizeDistribution::getIntegral()
{
    return this->distribution.sum();
}

Eigen::SparseVector<float, 0, int64_t> & InsertSizeDistribution::getDistributionVector()
{
    return this->distribution;
}

void InsertSizeDistribution::findMinProbability()
{
    // make sure that this is correct
    this->minValue = 1.0;
    
    for (uint32_t i = 0; i < this->distribution.nonZeros(); ++i)
        if (this->distribution.valuePtr()[i] > 0)
            this->minValue = std::min(this->distribution.valuePtr()[i], this->minValue);
    if (this->minValue == 1.0)
        this->minValue = 0.0;
}

float InsertSizeDistribution::getMinProbability()
{
    return this->minValue;
}
