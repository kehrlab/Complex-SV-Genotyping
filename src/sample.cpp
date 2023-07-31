#include "sample.hpp"
#include "bamFileHandler.hpp"
#include "custom_types.hpp"
#include "filter.hpp"
#include "genotypeResult.hpp"
#include "libraryDistribution.hpp"
#include "recordManager.hpp"
#include "variant.hpp"
#include "variantGenotyper.hpp"
#include <chrono>
#include <boost/filesystem.hpp>
#include <stdexcept>
#include <unordered_map>

Sample::Sample(
    Sample && s
): filename(s.filename), referenceFile(s.referenceFile), options(s.options), regionSampler(s.regionSampler), variants(s.variants)
{
    this->distributionDirectory = "";
    this->bamFileOpen = false;
}

Sample::Sample(
    std::string filename,
    SeqFileHandler & referenceFile, 
    ProgramOptions & options, 
    RegionSampler & regionSampler, 
    std::vector<complexVariant> & variants
) : filename(filename), referenceFile(referenceFile), options(options), regionSampler(regionSampler), variants(variants)
{
    this->distributionDirectory = "";
}

void Sample::open()
{
    openBamFile();
    createDistributionDirectory();
    calculateDefaultDistributions();
    closeBamFile();
}

void Sample::openBamFile()
{
    this->bamFile.open(this->filename, this->options);
    this->bamFileOpen = true;
}

void Sample::createDistributionDirectory()
{
    if (!this->options.isOptionOutputDistributions())
        return;

    this->distributionDirectory = this->filename + "_distributions";

    if (! boost::filesystem::create_directory(this->distributionDirectory))
    {
        if (! boost::filesystem::is_directory(this->distributionDirectory))
        {
            std::string msg = "Could not create directory '" + this->distributionDirectory + "', but flag -d was set. Make sure you have write permissions in the directory containing bam files.";
            std::runtime_error(msg.c_str());
        }
    }
}

void Sample::calculateDefaultDistributions()
{
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    std::unordered_map<std::string, TemplatePosition> insertPositions;
    if (this->options.isOptionUseWholeGenome()) { 
        insertPositions = this->bamFile.get_insert_size_positions();
    }
    else
    {
        std::vector<std::string> regionStrings;
        for (GenomicRegion & r : this->regionSampler.getSampledInsertRegions())
            regionStrings.push_back(r.getRegionString());
        insertPositions = this->bamFile.get_insert_size_positions(regionStrings);
    }
    createSampleDistribution(insertPositions);
    std::unordered_map<std::string, TemplatePosition>().swap(insertPositions);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    if (this->options.isOptionProfile())
        std::cout << this->filename << "\tTime taken to create default distributions: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms" << std::endl;
}

void Sample::createSampleDistribution(std::unordered_map<std::string, TemplatePosition> & insertPositions)
{
    if (this->options.isOptionGCCorrect() && !this->options.isOptionUseWholeGenome() && this->options.getRefFileName() != "")
    {
        this->sampleDistribution = LibraryDistribution(insertPositions, this->regionSampler.getSampledInsertRegions(), this->referenceFile);
    } else {
        this->sampleDistribution = LibraryDistribution(insertPositions);
    }
    if (this->options.isOptionOutputDistributions())
    {
        std::string defaultDistFileName = this->distributionDirectory + "/defaultDistribution.txt";
        this->sampleDistribution.writeDistribution(defaultDistFileName);
    }
}

LibraryDistribution & Sample::getLibraryDistribution()
{
    return this->sampleDistribution;
}

void Sample::closeBamFile()
{
    if (this->bamFileOpen)
    {
        this->bamFile.closeInputFile();
        this->bamFileOpen = false;
    }
}

int Sample::getFilterMargin()
{
    return ReadPairFilter::calculateBreakpointMargin(this->sampleDistribution.getInsertMean(), this->sampleDistribution.getInsertSD());
}

ContigInfo Sample::getContigInfo()
{
    return this->bamFile.getContigInfo();
}

std::unordered_map<std::string, int> Sample::getContigLengths()
{
    return this->bamFile.getContigLengths();
}

std::string Sample::getFileName()
{
    return this->filename;
}

void Sample::close()
{
    this->sampleDistribution = LibraryDistribution();
    closeBamFile();
}

