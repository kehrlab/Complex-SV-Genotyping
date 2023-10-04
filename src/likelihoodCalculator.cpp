#include "likelihoodCalculator.hpp"
#include "distributionConverter.hpp"
#include "genotypeResult.hpp"
#include "options.hpp"
#include "readTemplate.hpp"
#include "variant.hpp"
#include <chrono>
#include <limits>

LikelihoodCalculator::LikelihoodCalculator(RecordManager & recordManager, BamFileHandler & fileHandler, complexVariant & variant, LibraryDistribution & sampleDistribution, ProgramOptions & options):
    templates(recordManager.getTemplates()),
    variant(variant), 
    sampleDistribution(sampleDistribution),
    bamFileHandler(fileHandler),
    result(GenotypeResult(fileHandler.getFileName(), fileHandler.getContigInfo().cNames, options)),
    options(options),
    maxReadLength(recordManager.getMaxReadLength())
    {
        this->useInsertSizes = ! this->options.isOptionNoInsertSizes();
	    this->lowerInsertLimit = this->sampleDistribution.getInsertMean() - 2*this->sampleDistribution.getInsertSD();
	    this->upperInsertLimit = this->sampleDistribution.getInsertMean() + 2*this->sampleDistribution.getInsertSD();
        for (auto it = this->templates.begin(); it != this->templates.end(); ++it)
        {
            it->determineGroup(this->variant.getAllJunctions(), this->variant.getAllBreakpoints(), this->variant.getVariantRegions());
            it->determineLocationStrings();
        }
    };

void LikelihoodCalculator::calculateLikelihoods()
{
    this->result.initLikelihoods(this->genotypeNames);
    for (auto it = this->templates.begin(); it != this->templates.end(); ++it)
        if (it->isInterestingReadPair(this->filter)) 
            adjustLikelihoods(*it);
    this->templates = std::vector<ReadTemplate>();
    this->result.callGenotype();
}

void LikelihoodCalculator::adjustLikelihoods(ReadTemplate & readTemplate)
{
    float minValue = std::numeric_limits<float>::min();

    if (!readTemplate.isProperPair())
       return;

    int insertSize = readTemplate.getInsertSize();
    std::string orientation = readTemplate.getOrientation();
    bool split = readTemplate.containsSplitRead();
    bool spanning = readTemplate.containsSpanningRead();
    bool interChromosome = readTemplate.alignsAcrossChromosomes();
    std::string regionString = readTemplate.getRegionString();
    std::string junctionString = readTemplate.getJunctionString();
    std::string breakpointString = readTemplate.getBreakpointString();
    std::vector<int> mappingQualities = readTemplate.getMappingQualities();

    std::vector<std::vector<std::string>> rNamePairs = readTemplate.getInterChromosomeNames();

    if (spanning && this->options.isOptionNoSpanningReads())
	    return;
    if (split && this->options.isOptionNoSplitReads())
	    return;
    if (this->options.isOptionNoNormalReads())
	    if (!split && !spanning && orientation != "FF" && orientation != "RR" && insertSize > this->lowerInsertLimit && insertSize < this->upperInsertLimit)
		    return;

    std::vector<float> probabilities;
    bool insertSizeWithinLimits = false;
    for (auto dist : this->genotypeDistributions) {
	if (insertSize >= dist.getMinInsertSize() - 1000 && insertSize <= dist.getMaxInsertSize() + 1000)
		insertSizeWithinLimits = true;
        probabilities.push_back(std::max(dist.getProbability(insertSize, orientation, split, spanning, interChromosome, regionString, junctionString, breakpointString, rNamePairs, this->useInsertSizes), minValue));
    } 
    
    if (insertSizeWithinLimits)
    	this->result.storeEvidence(insertSize, orientation, split, spanning, interChromosome, regionString, junctionString, breakpointString, rNamePairs, mappingQualities);	
    
    this->result.addTemplateProbabilities(this->genotypeNames, probabilities, readTemplate.getTemplateWeight());
}

void LikelihoodCalculator::createInsertSizeDistributions()
{
    DistributionConverter distributionConverter(this->variant, this->sampleDistribution, this->maxReadLength, this->bamFileHandler, this->options);
    this->genotypeDistributions = distributionConverter.getGenotypeDistributions();
    this->genotypeNames = distributionConverter.getGenotypeNames();
    this->filter = distributionConverter.getReadPairFilter();
}

void LikelihoodCalculator::createInsertSizeDistributions(VariantProfile & variantProfile)
{
    DistributionConverter distributionConverter(variantProfile, this->sampleDistribution, this->bamFileHandler, this->options);
    this->genotypeDistributions = distributionConverter.getGenotypeDistributions();
    this->genotypeNames = distributionConverter.getGenotypeNames();
    this->filter = distributionConverter.getReadPairFilter();
}

GenotypeResult LikelihoodCalculator::getResult()
{
    return this->result;
}

std::vector<GenotypeDistribution> & LikelihoodCalculator::getDistributions()
{
    return this->genotypeDistributions;
}

std::vector<std::string> & LikelihoodCalculator::getGenotypeNames()
{
    return this->genotypeNames;
}

void LikelihoodCalculator::clearData()
{
    this->templates.clear();
    this->result = GenotypeResult();
}
