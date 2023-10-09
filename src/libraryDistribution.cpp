#include "libraryDistribution.hpp"
#include "custom_types.hpp"
#include "genomicRegion.hpp"
#include "readTemplate.hpp"
#include "seqan/sequence/sequence_shortcuts.h"
#include <cmath>
#include <fstream>
#include <unordered_map>

LibraryDistribution::LibraryDistribution()
{
    this->insertMean = 0;
    this->insertSD = 0;
    this->gcCorrection = false;
    this->sMin = 0;
    this->sMax = 1000;
    this->readLength = 0;
    this->numReadPairs = 0;
}

LibraryDistribution::LibraryDistribution(std::unordered_map<std::string, TemplatePosition> & positions)
{
    this->insertMean = 0;
    this->insertSD = 0;
    this->sMin = 0;
    this->sMax = 1000;
    this->readLength = 0;
    this->gcCorrection = false;
    this->numReadPairs = 0;

    this->distribution.resize(1001);
    for (int i = 0; i < this->distribution.size(); ++i)
        this->distribution[i] = std::vector<float>(100, 0);
    for (auto it : positions)
    {
        int s = it.second.end - it.second.begin + 1;
        if (s > 1000 || s < 0)
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

LibraryDistribution::LibraryDistribution(std::unordered_map<std::string, TemplatePosition> & insertPositions, std::vector<GenomicRegion> & regions, SeqFileHandler & refFileHandler)
{
    this->insertMean = 0;
    this->insertSD = 0;
    this->sMin = 0;
    this->sMax = 1000;
    this->gcCorrection = true;
    this->distribution.resize(1001);
    for (int i = 0; i < this->distribution.size(); ++i)
        this->distribution[i] = std::vector<float>(100, 0);
    this->uniformDistribution.resize(1001);
    for (int i = 0; i < this->uniformDistribution.size(); ++i)
        this->uniformDistribution[i] = std::vector<float>(100, 0);

    seqan::Dna5String gString = 'G';
    seqan::Dna5String cString = 'C';
    seqan::Dna5String base, seq;
    float gc, total, temp;
    int maxDist = 1000;

    std::vector<float> gcContents;

    // extract information from templates
    std::unordered_map<std::string, std::unordered_map<int, std::vector<int>>> templatePositions;
    for (auto it : insertPositions)
    {
        int s = it.second.end - it.second.begin + 1;
        if (s > 1000 || s < 0)
            continue;
        this->insertSizes.push_back(s);
        ++this->numReadPairs;

        if (templatePositions.find(it.second.chr) == templatePositions.end())
        {
            std::unordered_map<int, std::vector<int>> tempPositions;
            std::vector<int> endPos;
            endPos.push_back(s);
            tempPositions[it.second.begin] = endPos;
            templatePositions[it.second.chr] = tempPositions;
        } 
        else 
        {
            if (templatePositions[it.second.chr].find(it.second.begin) == templatePositions[it.second.chr].end())
            {
                std::vector<int> endPos;
                endPos.push_back(s);
                templatePositions[it.second.chr][it.second.begin] = endPos;
            } 
            else 
            {
                templatePositions[it.second.chr][it.second.begin].push_back(s);
            }
        }
    }

    // traverse all regions
    for (GenomicRegion & region : regions)
    {
        if (!region.sequenceInMemory())
            region.readSequence(refFileHandler);
        seq = region.getSequence();

        std::unordered_map<int, std::vector<int>> positions;
        if (templatePositions.find(region.getReferenceName()) != templatePositions.end())
            positions = templatePositions[region.getReferenceName()];
        
        for (int i = 0; i < seqan::length(seq); ++i)
        {
            if (i == 0)
            {
                total = 0;
                gc = 0;
                for (int s = 0; s <= maxDist; ++s)
                {
                    if ((i + s) >= seqan::length(seq))
                        break;
                    base = seq[i+s];
                    ++total;
                    if (base == gString || base == cString)
                        ++gc;
                    gcContents.push_back(gc/total);
                }
            } else {
                total = 1;
                base = seq[i - 1];
                for (int j = 0; j < gcContents.size() - 1; ++j)
                {
                    if ((i + j) >= (seqan::length(seq) - 1))
                        break;
                    ++total;
                    temp = gcContents[j + 1] * total;
                    if (base == gString || base == cString)
                        temp -= 1;
                    gcContents[j] = temp / (total - 1);
                }
                if ((i + maxDist) < seqan::length(seq))
                {
                    total = maxDist + 1;
                    temp = gcContents[gcContents.size() - 1] * total;
                    if (base == gString || base == cString)
                            temp -= 1;
                    base = seq[i + maxDist];
                    if (base == gString || base == cString)
                            temp += 1;
                    gcContents[gcContents.size() - 1] = temp / total;
                }
            }

            for (int j = 0; j < gcContents.size(); ++j)
            {
                if ((i+j) >= seqan::length(seq))
                    break;
                this->uniformDistribution[j][calculateGCIndex(gcContents[j])] += 1;
            }
            if (positions.find(region.getRegionStart() + i) != positions.end())
            {
                for (int s : positions[region.getRegionStart() + i])	
			this->distribution[s][calculateGCIndex(gcContents[s])] += 1;
            }
        }
        gcContents.erase(gcContents.begin(), gcContents.end());
        region.clearSequence();
    }
    smoothDistribution(this->distribution, 5); 
    scaleDistribution(this->distribution);
    
    smoothDistribution(this->uniformDistribution, 5);
    scaleDistribution(this->uniformDistribution);
    
    createMarginalDistributions();
    calculateCorrectionFactors();

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
    for (int y = 0; y < distribution[0].size(); ++y)
    {
        smoothed_values.resize(distribution.size());
        for (int x = 0; x < distribution.size(); ++x)
        {
            float smoothed = 0;
            for (int i = 0; i < kernel.size(); ++i) 
            {
                if (x+i-(kernelSize-1)/2 < 0 || x+i-(kernelSize-1)/2 >= distribution.size())
                    continue;
                smoothed += kernel[i] * distribution[x+i-(kernelSize-1)/2][y];
            }
            smoothed_values[x] = smoothed;
        }
        for (int x = 0; x < distribution.size(); ++x)
            distribution[x][y] = smoothed_values[x];
    }

    if (!this->gcCorrection)
        return;

    // smooth in gc-direction
    for (int x = 0; x < distribution.size(); ++x)
    {
        smoothed_values.resize(distribution[0].size());
        for (int y = 0; y < distribution[x].size(); ++y)
        {
            float smoothed = 0;
            for (int i = 0; i < kernel.size(); ++i) 
            {
                if (y+i-(kernelSize-1)/2 < 0 || y+i-(kernelSize-1)/2 >= distribution[x].size())
                    continue;
                smoothed += kernel[i] * distribution[x][y+i-(kernelSize-1)/2];
            }
            smoothed_values[y] = smoothed;
        }
	for (int y = 0; y < distribution[x].size(); ++y)
        	distribution[x][y] = smoothed_values[y];
    }
}

void LibraryDistribution::scaleDistribution(std::vector<std::vector<float>> & distribution)
{
    double pTotal = 0.0;
    for (int x = 0; x < distribution.size(); ++x)
        for (int y = 0; y < distribution[x].size(); ++y)
            pTotal += distribution[x][y];
    for (int x = 0; x < distribution.size(); ++x)
        for (int y = 0; y < distribution[x].size(); ++y)
            distribution[x][y] = (float) (distribution[x][y] / pTotal);
}

void LibraryDistribution::createMarginalDistributions()
{
    this->insertDistribution.resize(this->distribution.size());
    this->gcDistribution.resize(this->distribution[0].size());

    // insert size distribution
    for (int i = 0; i < this->distribution.size(); ++i)
    {
        this->insertDistribution[i] = 0.0;
        for (int j = 0; j < this->distribution[i].size(); ++j)
            this->insertDistribution[i] += this->distribution[i][j];
    }
    calculateInsertStats();

    // gc distribution
    for (int i = 0; i < this->distribution[0].size(); ++i)
    {
        this->gcDistribution[i] = 0.0;
        for (int j = 0; j < this->distribution.size(); ++j)
            this->gcDistribution[i] += this->distribution[j][i];
    }
}

void LibraryDistribution::calculateCorrectionFactors()
{
    this->correctionFactors.resize(1001);
    for (int i = 0; i < this->correctionFactors.size(); ++i)
	    this->correctionFactors[i] = std::vector<float>(100, 0);
    
    for (int s = 0; s < this->distribution.size(); ++s)
    {
        float pMargin = 0;
        for (int j = 0; j < this->distribution[s].size(); ++j)
            pMargin += this->uniformDistribution[s][j];
        
        // calculate
	if (this->insertDistribution[s] == 0)
		continue;	
        float marginalFactor = pMargin / this->insertDistribution[s];
        for (int j = 0; j < this->correctionFactors[s].size(); ++j)
	{
		if (this->uniformDistribution[s][j] == 0)
			continue;
            this->correctionFactors[s][j] = (this->distribution[s][j] / this->uniformDistribution[s][j]) * marginalFactor;
	}
    }
}

float LibraryDistribution::getProbability(int s, float gc)
{
    if (s >= 0 && s <= 1000 && gc >= 0 && gc <= 1)
        return this->distribution[s][calculateGCIndex(gc)];
    return 0.0;
}

float LibraryDistribution::getProbability(int s)
{
    if (s >= 0 && s <= 1000)
        return this->insertDistribution[s];
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
        for (int s = 0; s < this->distribution.size(); ++s) {
            if (this->gcCorrection) {
		f << s << "\t" << 0.01 << "\t" << this->distribution[s][0] << "\t" << this->correctionFactors[s][0] << std::endl;
                for (int j = 1; j < this->distribution[s].size(); ++j)
                    f << s << "\t" << (j+1) * 0.01 << "\t" << this->distribution[s][j] << "\t" << this->correctionFactors[s][j] << std::endl;
            } else {
		    f << s << "\t" << this->distribution[s][0] << std::endl;
	    }
        }
    } else {
        std::cerr << "Could not open file: " << filename << "for writing." << std::endl;
    }
    f.close();
}

float LibraryDistribution::getCorrectionFactor(int s, float gc)
{	
    return this->correctionFactors[s][calculateGCIndex(gc)];
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