#include "libraryDistribution.hpp"

LibraryDistribution::LibraryDistribution()
{
    this->insertMean = 0;
    this->insertSD = 0;
    this->gcCorrection = false;
    this->sMin = 1;
    this->sMax = 1000;
    this->readLength = 0;
    this->numReadPairs = 0;
}

LibraryDistribution::LibraryDistribution(std::unordered_map<std::string, TemplatePosition> & positions)
{
    this->insertMean = 0;
    this->insertSD = 0;
    this->sMin = 1;
    this->sMax = 1000;
    this->readLength = 0;
    this->gcCorrection = false;
    this->numReadPairs = 0;

    this->distribution.resize(this->sMax - this->sMin + 1);
    for (uint32_t i = 0; i < this->distribution.size(); ++i)
        this->distribution[i] = std::vector<float>(100, 0);
    for (auto it : positions)
    {
        int s = it.second.end - it.second.begin + 1;
        if (s > this->sMax || s < this->sMin)
            continue;
        ++this->numReadPairs;
        this->insertSizes.push_back(s);
        this->distribution[s][0] += 1;
    }
    smoothDistribution(this->distribution, 10);
    smoothDistribution(this->distribution, 20);
    smoothDistribution(this->distribution, 5);
    scaleDistribution(this->distribution);
    createMarginalDistributions();
    std::vector<int>().swap(this->insertSizes);
}

LibraryDistribution::LibraryDistribution(int sMin, int sMax, int readLength, std::vector<double> & values)
{
    this->insertMean = 0;
    this->insertSD = 0;
    this->gcCorrection = false;
    this->sMin = sMin;
    this->sMax = sMax;
    this->readLength = readLength;
    this->numReadPairs = 0;

    this->distribution.resize(this->sMax - this->sMin + 1);
    for (uint32_t i = 0; i < this->distribution.size(); ++i)
        this->distribution[i] = std::vector<float>(100, 0);
    for (uint32_t i = 0; i < values.size(); ++i)
    {
        uint32_t nR = (int) values[i];
        this->numReadPairs += (int) values[i];
        this->distribution[i][0] += (float) values[i];
        for (uint32_t j = 0; j < nR; ++j)
            this->insertSizes.push_back(this->sMin + i);
    }
    smoothDistribution(this->distribution, 10);
    smoothDistribution(this->distribution, 20);
    smoothDistribution(this->distribution, 5);
    scaleDistribution(this->distribution);
    createMarginalDistributions();
    std::vector<int>().swap(this->insertSizes);
}

int LibraryDistribution::calculateGCIndex(float gcContent)
{
    return std::min((int) (gcContent / 0.01), 99);
}

void LibraryDistribution::smoothDistribution(std::vector<std::vector<float>> & distribution, int kernelSize)
{
    std::vector<float> kernel(kernelSize, 0);
    float sigma = (float) kernelSize / 5;
    for (int i = - (kernelSize-1)/2; i <= (kernelSize-1)/2; ++i)
        kernel[i+(kernelSize-1)/2] = (1 / (std::sqrt(2 * M_PI) * sigma)) * std::exp(- (i*i)/(2*sigma*sigma));

    
    // smooth first in s-direction
    std::vector<float> smoothed_values;
    for (uint32_t y = 0; y < distribution[0].size(); ++y)
    {
        smoothed_values.resize(distribution.size());
        for (int x = 0; x < (int) distribution.size(); ++x)
        {
            float smoothed = 0;
            for (int i = 0; i < (int) kernel.size(); ++i) 
            {
                if (x + i - (kernelSize - 1)/2 < 0 || x + i - (kernelSize - 1)/2 >= (int) distribution.size())
                    continue;
                smoothed += kernel[i] * distribution[x+i-(kernelSize-1)/2][y];
            }
            smoothed_values[x] = smoothed;
        }
        for (uint32_t x = 0; x < distribution.size(); ++x)
            distribution[x][y] = smoothed_values[x];
    }

    if (!this->gcCorrection)
        return;

    // smooth in gc-direction
    for (uint32_t x = 0; x < distribution.size(); ++x)
    {
        smoothed_values.resize(distribution[0].size());
        for (int y = 0; y < (int) distribution[x].size(); ++y)
        {
            float smoothed = 0;
            for (int i = 0; i < (int) kernel.size(); ++i) 
            {
                if (y+i-(kernelSize-1)/2 < 0 || y+i-(kernelSize-1)/2 >= (int) distribution[x].size())
                    continue;
                smoothed += kernel[i] * distribution[x][y+i-(kernelSize-1)/2];
            }
            smoothed_values[y] = smoothed;
        }
	    for (uint32_t y = 0; y < distribution[x].size(); ++y)
        	distribution[x][y] = smoothed_values[y];
    }
}

void LibraryDistribution::scaleDistribution(std::vector<std::vector<float>> & distribution)
{
    double pTotal = 0.0;
    for (uint32_t x = 0; x < distribution.size(); ++x)
        for (uint32_t y = 0; y < distribution[x].size(); ++y)
            pTotal += distribution[x][y];
    for (uint32_t x = 0; x < distribution.size(); ++x)
        for (uint32_t y = 0; y < distribution[x].size(); ++y)
            distribution[x][y] = (float) (distribution[x][y] / pTotal);
}

void LibraryDistribution::createMarginalDistributions()
{
    this->insertDistribution.resize(this->distribution.size());
    this->gcDistribution.resize(this->distribution[0].size());

    // insert size distribution
    for (uint32_t i = 0; i < this->distribution.size(); ++i)
    {
        this->insertDistribution[i] = 0.0;
        for (uint32_t j = 0; j < this->distribution[i].size(); ++j)
            this->insertDistribution[i] += this->distribution[i][j];
    }
    calculateInsertStats();

    // gc distribution
    for (uint32_t i = 0; i < this->distribution[0].size(); ++i)
    {
        this->gcDistribution[i] = 0.0;
        for (uint32_t j = 0; j < this->distribution.size(); ++j)
            this->gcDistribution[i] += this->distribution[j][i];
    }
}

void LibraryDistribution::calculateCorrectionFactors()
{
    this->correctionFactors.resize(this->sMax - this->sMin + 1);
    for (uint32_t i = 0; i < this->correctionFactors.size(); ++i)
	    this->correctionFactors[i] = std::vector<float>(100, 0);
    
    for (uint32_t s = 0; s < this->distribution.size(); ++s)
    {
        float pMargin = 0;
        for (uint32_t j = 0; j < this->distribution[s].size(); ++j)
            pMargin += this->uniformDistribution[s][j];
        
        // calculate
        if (this->insertDistribution[s] == 0)
            continue;	
        float marginalFactor = pMargin / this->insertDistribution[s];
        for (uint32_t j = 0; j < this->correctionFactors[s].size(); ++j)
        {
            if (this->uniformDistribution[s][j] == 0)
                continue;
            this->correctionFactors[s][j] = (this->distribution[s][j] / this->uniformDistribution[s][j]) * marginalFactor;
        }
    }
}

float LibraryDistribution::getProbability(int s, float gc)
{
    if (s >= this->sMin && s <= this->sMax && gc >= 0 && gc <= 1)
        return this->distribution[s-this->sMin][calculateGCIndex(gc)];
    return 0.0;
}

float LibraryDistribution::getProbability(int s)
{
    if (s >= this->sMin && s <= this->sMax)
        return this->insertDistribution[s - this->sMin];
    return 0.0;
}

float LibraryDistribution::getProbability(float gc)
{
    if (gc >= 0 && gc <= 1)
        return this->gcDistribution[calculateGCIndex(gc)];
    return 0.0;
}

void LibraryDistribution::writeDistribution(std::string filename)
{
    std::ofstream f(filename);
    if (f.is_open())
    {
	    if (this->gcCorrection) {
        	f << "InsertSize\tGC\tProbability\tCorrection" << std::endl;
	    } else {
		    f << "InsertSize\tProbability" << std::endl;
	    }
        for (int s = 0; s < (int) this->insertDistribution.size(); ++s) {
            if (this->gcCorrection) {
		        f << s + this->sMin << "\t" << 0.01 << "\t" << this->distribution[s][0] << "\t" << this->correctionFactors[s][0] << std::endl;
                for (uint32_t j = 1; j < this->distribution[s].size(); ++j)
                    f << s + this->sMin << "\t" << (j+1) * 0.01 << "\t" << this->distribution[s][j] << "\t" << this->correctionFactors[s][j] << std::endl;
            } else {
		        f << s << "\t" << this->insertDistribution[s] << std::endl;
	        }
        }
    } else {
        std::cerr << "Could not open file: " << filename << " for writing." << std::endl;
    }
    f.close();
}

float LibraryDistribution::getCorrectionFactor(int s, float gc)
{	
    return this->correctionFactors[s - this->sMin][calculateGCIndex(gc)];
}

void LibraryDistribution::calculateInsertStats()
{
    this->insertMean = 0;
    this->insertSD = 0;

    for (auto it : this->insertSizes)
        this->insertMean += it;
    this->insertMean /= this->insertSizes.size();

    for (auto it : this->insertSizes)
        this->insertSD += (it - this->insertMean) * (it-this->insertMean);
    this->insertSD = std::sqrt(this->insertSD / this->insertSizes.size());
}

float & LibraryDistribution::getInsertMean()
{
    return this->insertMean;
}

float & LibraryDistribution::getInsertSD()
{
    return this->insertSD;
}

bool LibraryDistribution::gcCorrectionPossible()
{
    return this->gcCorrection;
}

int LibraryDistribution::drawInsertSize()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<int> distrib(this->insertDistribution.begin(), this->insertDistribution.end());
    return distrib(gen);
}

int & LibraryDistribution::getMinInsert()
{
    return this->sMin;
}

int & LibraryDistribution::getMaxInsert()
{
    return this->sMax;
}

std::vector<float> & LibraryDistribution::getInsertDistribution()
{
    return this->insertDistribution;
}

int & LibraryDistribution::getNumReadPairs()
{
    return this->numReadPairs;
}

int & LibraryDistribution::getReadLength()
{
    return this->readLength;
}