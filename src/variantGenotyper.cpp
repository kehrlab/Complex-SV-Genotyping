#include "variantGenotyper.hpp"
#include "bamFileHandler.hpp"
#include "libraryDistribution.hpp"
#include "likelihoodCalculator.hpp"

VariantGenotyper::VariantGenotyper(
    complexVariant & variant,
    Sample & sample,
    ProgramOptions & options
) : variant(variant), bamFileName(bamFileName), sampleName(sample.getSampleName()), sampleDistribution(sample.getLibraryDistribution()), options(options)
{}

VariantGenotyper::VariantGenotyper(
    complexVariant & variant, 
    std::string bamFileName,
    std::string sampleName,
    LibraryDistribution & sampleDistribution,
    ProgramOptions & options
) : variant(variant), bamFileName(bamFileName), sampleName(sampleName), sampleDistribution(sampleDistribution), options(options)
{}

void VariantGenotyper::genotype()
{
    BamFileHandler bamFileHandler(this->bamFileName, this->options);
    bamFileHandler.setRegions(variant.calculateAssociatedRegions(this->sampleDistribution));
    this->bamRecords.setReadPairs(bamFileHandler.get_read_pairs());

    std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> chromosomeStructures;
    for (auto & allele : this->variant.getAlleles())
        chromosomeStructures[allele.getName()] = allele.getChromosomeStructures(this->variant.getAllBreakpoints());

    LikelihoodCalculator likelihoodCalculator(bamRecords, bamFileHandler, variant, this->sampleName, this->sampleDistribution, options);
    likelihoodCalculator.createInsertSizeDistributions(chromosomeStructures);
    this->bamRecords.reset();
    likelihoodCalculator.calculateLikelihoods(chromosomeStructures);
    this->variantGenotype = likelihoodCalculator.getResult();
    if (this->options.isOptionOutputDistributions())
        writeInsertSizeDistributions(likelihoodCalculator);
    this->variantGenotype.clearData();
    bamFileHandler.closeInputFile();
}

void VariantGenotyper::genotype(VariantProfile & variantProfile)
{
    BamFileHandler bamFileHandler(this->bamFileName, this->options);
    bamFileHandler.setRegions(variant.calculateAssociatedRegions(this->sampleDistribution));
    this->bamRecords.setReadPairs(bamFileHandler.get_read_pairs());
    LikelihoodCalculator likelihoodCalculator(bamRecords, bamFileHandler, variantProfile.getVariant(), this->sampleName, this->sampleDistribution, options);
    likelihoodCalculator.createInsertSizeDistributions(variantProfile);
    this->bamRecords.reset();
    likelihoodCalculator.calculateLikelihoods(variantProfile);
    this->variantGenotype = likelihoodCalculator.getResult();
    if (this->options.isOptionOutputDistributions())
        writeInsertSizeDistributions(likelihoodCalculator);
    this->variantGenotype.clearData();
    bamFileHandler.closeInputFile();
}

void VariantGenotyper::writeInsertSizeDistributions(LikelihoodCalculator & likelihoodCalculator)
{
    if (!this->options.isOptionOutputDistributions())
        return;

    std::string variantName = "unnamed";
    if (this->variant.getName() != "")
        variantName = this->variant.getName();
    
    std::string directory = this->bamFileName + "_distributions/" + variantName;
    if (! boost::filesystem::create_directory(directory)) {
        if (! boost::filesystem::is_directory(directory)) {
            std::cerr << "Could not create variant distribution directory" << std::endl;
            return;
        }
    }
    directory += "/";

    
    // write theoretical distributions
    std::vector<std::string> & genotypeNames = likelihoodCalculator.getGenotypeNames();
    std::vector<GenotypeDistribution> & distributions = likelihoodCalculator.getDistributions();

    for (int i = 0; i < genotypeNames.size(); ++i)
    {
        std::string prefix = directory + "numericalDistribution";
        
        std::string gtName = genotypeNames[i];
        int slashIndex = gtName.find("/");
        gtName.replace(slashIndex, 1, "-");

        if (variantName != "")
            prefix.append("_" + variantName);
        prefix.append("_" + gtName);
        distributions[i].writeDistributionBinned(prefix);
    }

    // write observed distribution
    std::string observedDistributionPrefix = directory + "observedDistribution_" + variantName;
    this->variantGenotype.writeEvidence(observedDistributionPrefix);

    // write bootstrap result
    std::string bootstrapPrefix = directory + variantName;
    this->variantGenotype.writeBootstrapData(bootstrapPrefix);
}

GenotypeResult VariantGenotyper::getResult()
{
    return this->variantGenotype;
}

std::string VariantGenotyper::getVariantName()
{
    return this->variant.getName();
}

std::string VariantGenotyper::getSampleName()
{
    return this->bamFileName;
}