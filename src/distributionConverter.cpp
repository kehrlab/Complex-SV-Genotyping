#include "distributionConverter.hpp"
#include "allele.hpp"
#include "custom_types.hpp"
#include "genomicRegion.hpp"
#include "genotypeDistribution.hpp"
#include "libraryDistribution.hpp"
#include "options.hpp"
#include "readTemplate.hpp"
#include "record.hpp"
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <fstream>

DistributionConverter::DistributionConverter()
{
    this->readLength = 0;
}

DistributionConverter::DistributionConverter(complexVariant variant, LibraryDistribution sampleDistribution, int readLength, BamFileHandler & fileHandler, ProgramOptions & options)
{
    this->readLength = readLength;
    this->sampleDistribution = sampleDistribution;
    this->variant = variant;
    this->filter = ReadPairFilter(this->variant.getAllBreakpoints(), 500, 0);//this->sampleDistribution.getInsertMean(), this->sampleDistribution.getInsertSD());
    this->options = options;

    this->filename = fileHandler.getFileName();

    // this->lowerInsertLimit = this->sampleDistribution.getInsertMean() - 2*this->sampleDistribution.getInsertSD();
    // this->upperInsertLimit = this->sampleDistribution.getInsertMean() + 2*this->sampleDistribution.getInsertSD();
    this->lowerInsertLimit = 1;
    this->upperInsertLimit = 1000;

    // createVariantMaps();
    initDistributions(fileHandler, options.getDistributionMode());
    createDistributions();
    determineDifficulty();
}

DistributionConverter::DistributionConverter(VariantProfile & variantProfile, LibraryDistribution & sampleDistribution, BamFileHandler & fileHandler, ProgramOptions & options)
{
    this->filter = variantProfile.getFilter();
    this->options = options;
    this->filename = fileHandler.getFileName();
    this->filter = variantProfile.getFilter();

    std::vector<std::string> cNames = fileHandler.getContigInfo().cNames;
    auto tempDists = variantProfile.calculateGenotypeDistributions(sampleDistribution, 0.001);
    for (auto & dist : tempDists)
    {
        dist.second.setPossibleContigs(cNames);
        this->genotypeDistributions.push_back(dist.second);
        this->genotypeNames.push_back(dist.first);
    }
    determineDifficulty();
}

void DistributionConverter::initDistributions(BamFileHandler & fileHandler, int distributionMode)
{
    std::vector<Allele> & alleles = this->variant.getAlleles();
    std::vector<std::string> cNames = fileHandler.getContigInfo().cNames;

    for (auto it = alleles.begin(); it != alleles.end(); ++it)
    {
        this->alleleNames.push_back(it->getName());
        GenotypeDistribution tempDistribution(cNames, distributionMode);
        this->alleleDistributions.push_back(tempDistribution);
    }
}


void DistributionConverter::createDistributions()
{
    createVariantDistributions();
    createMixedDistributions();
    scaleDistributions();
}

void DistributionConverter::determineDifficulty()
{
    if (this->options.getEstimateDiffQual() == 0)
        return;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Variant: " << this->variant.getName() << ";\tFile: " << this->filename << std::endl;
    std::cout << "Expected number of reads required to reach quality of " << this->options.getEstimateDiffQual() << ": " << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Ground Truth\t#Reads" << std::endl;
    for (int i = 0; i < this->genotypeDistributions.size(); ++i) 
    {
        std::cout << this->genotypeNames[i] << "\t";
        int requiredReads = 0;
        std::string altName;
        for (int j = 0; j < this->genotypeDistributions.size(); ++j)
        {
            if (j != i)
            {
                int tempReads = (int) (this->options.getEstimateDiffQual() / this->genotypeDistributions[i].calculateKLD(this->genotypeDistributions[j]));
                if (tempReads > requiredReads)
                {
                    altName = this->genotypeNames[j];
                    requiredReads = tempReads;
                }
            }
        }
        std::cout << requiredReads << " (" << altName << ")" << std::endl;
        
    }
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    return;
}

void DistributionConverter::createVariantDistributions()
{
    for (int i = 0; i < this->variant.getAlleles().size(); ++i) {
	    createAlleleDistribution(this->alleleDistributions[i], this->variant.getAlleles()[i]);
        if (this->options.isOptionOutputDistributions())
        {
            std::string alleleFileName = this->filename + "_distributions/" + this->variant.getName() +"/" + this->variant.getAlleles()[i].getName();
            this->alleleDistributions[i].writeDistributionBinned(alleleFileName);
        }
    }
}

void DistributionConverter::createAlleleDistribution(GenotypeDistribution & distribution, Allele & allele)
{
    for (std::string cName : allele.getChromosomeNames()) 
        addVariantRegions(cName, distribution, allele, this->sampleDistribution.gcCorrectionPossible());
    if (this->options.isOptionSample() && !this->sampleDistribution.gcCorrectionPossible())
        distribution.smoothDistribution();
}

void DistributionConverter::addVariantRegions(std::string cName, GenotypeDistribution & distribution, Allele & allele, bool correctionEnabled)
{
    if (correctionEnabled)
        addVariantRegionsWithGCCorrection(cName, distribution, allele);
    else
        addVariantRegionsWithoutGCCorrection(cName, distribution, allele);
}

void DistributionConverter::addVariantRegionsWithoutGCCorrection(std::string cName, GenotypeDistribution & distribution, Allele & allele)
{
    VariantMapManager & chrMaps = allele.getChromosomeMap(cName);
    std::vector<VariantMap> maps = chrMaps.getMaps();
    int maxInsertSize = 1000;
    int minInsertSize = 0;
    float insertProbability = 0;
    int status = 0;

    // init random generator
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int j = 0; j < maps.size(); ++j)
    {
	    VariantMap map = maps[j];
        int mapLength = map.totalLength;

        if (this->options.isOptionSample())
        {
            int NR = (int) (this->options.getCoverage() * mapLength / (2 * this->readLength)) + 1;
            int iterations = 0;
            // std::cout << "Sampling " << NR << " reads to achieve a coverage of " << this->options.getCoverage() << std::endl;
            int sampledReads = 0;
            int position = 0;
            int insertSize = 0;
            
            std::uniform_int_distribution<> positionDistribution(0, mapLength - 1);

            while (sampledReads < NR) 
            {
                ++iterations;
                position = positionDistribution(gen);
                insertSize = this->sampleDistribution.drawInsertSize();
                if (position >= mapLength || position + insertSize >= mapLength)
                    continue;
                addSimulatedTemplateToDistribution(map, position, insertSize, 1.0, distribution, allele);
                ++sampledReads;
            }
        } else {
            for (unsigned i = 0; i < mapLength; ++i)
            {
                for (int s = maxInsertSize; s >= minInsertSize; --s)
                {
                    if (s >= mapLength || i+s >= mapLength)
                        continue;
                    insertProbability = this->sampleDistribution.getProbability(s);
                    if (insertProbability == 0)
                        continue;

                    status = addSimulatedTemplateToDistribution(map, i, s, insertProbability, distribution, allele);
                    if (status < 0)
                        break;
                }
            }
        }
    }
}

void DistributionConverter::addVariantRegionsWithGCCorrection(std::string cName, GenotypeDistribution & distribution, Allele & allele)
{
    VariantMapManager & chrMaps = allele.getChromosomeMap(cName);
    std::vector<VariantMap> maps = chrMaps.getMaps();
    int maxInsertSize = 1000;
    int minInsertSize = 0;

    for (int j = 0; j < maps.size(); ++j)
    {
	    VariantMap map = maps[j]; // the map contains the allele sequence
        int mapLength = map.totalLength;
        std::vector<float> gcContentVector;

        for (unsigned i = 0; i < mapLength; ++i)
        {
            map.getGCContentVector(gcContentVector, i, minInsertSize, maxInsertSize);
            for (int s = maxInsertSize; s >= minInsertSize; --s)
            {
                if (s >= mapLength || i+s >= mapLength)
                    continue;
                float insertProbability = sampleDistribution.getProbability(s);
                if (insertProbability == 0)
                    continue;

                ReadTemplate sT = map.simulateTemplate(i, s, this->readLength);
                sT.findSpanningReads(this->variant.getAllBreakpoints());
                if (sT.containsSuspectedSplit()) 
			        sT.findSplitReads(allele.getNovelJunctions());
                sT.determineOverlappingRegions(this->variant.getVariantRegions());
                sT.determineLocationStrings();

                std::string orientation = sT.getOrientation();
                int insertSize = sT.getInsertSize();
                
                bool split = sT.containsSplitRead();
                if (split && this->options.isOptionNoSplitReads()) {
                    continue;
		        }

                bool spanning = sT.containsSpanningRead();
                if (spanning && this->options.isOptionNoSpanningReads()) {
                    continue;
		        }

                if (this->options.isOptionNoNormalReads())
                    if (!split && !spanning && orientation != "FF" && orientation != "RR" && insertSize > this->lowerInsertLimit && insertSize < this->upperInsertLimit) { 
                        continue;
		        }

                bool interChromosome = sT.alignsAcrossChromosomes();
                std::string regionString = sT.getRegionString();
                std::string junctionString = sT.getJunctionString();
                std::string breakpointString = sT.getBreakpointString();

                int idx = s - minInsertSize;
                float gcCorrectionFactor = extractGCCorrectionFactor(gcContentVector[idx], s);
                insertProbability *= gcCorrectionFactor;

                if (insertProbability == 0) {
                    continue;
		}

                std::vector<std::vector<std::string>> rNamePairs = sT.getInterChromosomeNames();

                if (sT.isInterestingReadPair(this->filter)) {
                    distribution.addInsertSizeProbability(insertSize, orientation, split, spanning, interChromosome, regionString, junctionString, breakpointString, rNamePairs, insertProbability);
             
	       	}else{
                    break;
                }
            }
        }
    }
}

int DistributionConverter::addSimulatedTemplateToDistribution(VariantMap & map, int i, int s, float insertProbability, GenotypeDistribution & distribution, Allele & allele)
{
    ReadTemplate sT = map.simulateTemplate(i, s, this->readLength);
    sT.findSpanningReads(this->variant.getAllBreakpoints());
    if (sT.containsSuspectedSplit()) 
        sT.findSplitReads(allele.getNovelJunctions());
    sT.determineOverlappingRegions(this->variant.getVariantRegions());
    sT.determineLocationStrings();

    std::string orientation = sT.getOrientation();
    int insertSize = sT.getInsertSize();

    bool split = sT.containsSplitRead();
    bool spanning = sT.containsSpanningRead();
    bool interChromosome = sT.alignsAcrossChromosomes();
    std::string regionString = sT.getRegionString();
    std::string junctionString = sT.getJunctionString();
    std::string breakpointString = sT.getBreakpointString();

    std::vector<std::vector<std::string>> rNamePairs = sT.getInterChromosomeNames();

    if (split && this->options.isOptionNoSplitReads())
        return 0;
    if (spanning && this->options.isOptionNoSpanningReads())
        return 0;
    if (this->options.isOptionNoNormalReads())
        if (!split && !spanning && orientation != "FF" && orientation != "RR" && insertSize > this->lowerInsertLimit && insertSize < this->upperInsertLimit)
            return 0;
    if (insertProbability == 0)
        return 0;
    
    if (sT.isInterestingReadPair(this->filter)) {
        distribution.addInsertSizeProbability(insertSize, orientation, split, spanning, interChromosome, regionString, junctionString, breakpointString, rNamePairs, insertProbability);
        return 1;
    } else {
        return -1;
    }
    return 0;
}

float DistributionConverter::extractGCCorrectionFactor(float gcContent, int insertSize)
{
    float gcCorrectionFactor = 1.0;
    if (gcContent == 0)
        return gcCorrectionFactor;
    return this->sampleDistribution.getCorrectionFactor(insertSize, gcContent);
}

void DistributionConverter::createMixedDistributions()
{
    double eps = 0.001;
    for (int i = 0; i < this->alleleNames.size(); ++i)
    {
        GenotypeDistribution tempDistribution;
        for (int j = 0; j < this->alleleNames.size(); ++j)
        {
            if (j == i)
                continue;
            tempDistribution += this->alleleDistributions[j];
        }
        this->genotypeNames.push_back(alleleNames[i] + "/" + alleleNames[i]);
        this->genotypeDistributions.push_back(
            (1.0-eps) * alleleDistributions[i] + eps * tempDistribution
        );

        for (int j = 0; j < i; ++j)
        {
            this->genotypeNames.push_back(alleleNames[j] + "/" + alleleNames[i]);
            GenotypeDistribution mixedDist = this->alleleDistributions[j] + this->alleleDistributions[i];
            if (this->alleleNames.size() == 2)
            {
                this->genotypeDistributions.push_back(mixedDist);
            }
            else
            {
                GenotypeDistribution otherDist;
                for (int k = 0; k < this->alleleNames.size(); ++k)
                {
                    if (k == j || k == i)
                        continue;
                    otherDist += this->alleleDistributions[k];
                }
                this->genotypeDistributions.push_back((1-eps) / 2 * mixedDist + (eps /(this->alleleNames.size() - 2)) * otherDist);
            }
        }
    }
}

void DistributionConverter::scaleDistributions()
{
    for (int i = 0; i < this->genotypeDistributions.size(); ++i)
        this->genotypeDistributions[i].scaleDistribution();
    setMinProbabilities(determineMinProbability());
}

float DistributionConverter::determineMinProbability()
{
    float minProbability = 1;
    for (auto & dist: this->genotypeDistributions)
        if (dist.getMinProbability() > 0 && dist.getMinProbability() < minProbability)
            minProbability = dist.getMinProbability();
    return minProbability;
}

void DistributionConverter::setMinProbabilities(float minProbability)
{
    for (auto & dist: this->genotypeDistributions)
        dist.setMinProbability(minProbability);
}

ReadPairFilter & DistributionConverter::getReadPairFilter()
{
    return this->filter;
}

std::vector<GenotypeDistribution> & DistributionConverter::getGenotypeDistributions()
{
    return this->genotypeDistributions;
}

std::vector<std::string> & DistributionConverter::getGenotypeNames()
{
    return this->genotypeNames;
}
