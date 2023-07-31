#include "options.hpp"
#include "genomicRegion.hpp"
#include <fstream>
#include <stdexcept>
#include <iostream>

ProgramOptions::ProgramOptions()
{
    this->outFile = "";
    this->variantFile = "";
    this->refFile = "";
    this->output = false;
    this->verbose = false;
    this->wholeGenome = false;
    this->outputDistributions = false;
    this->estimateDiffQ = 0;
    this->nThreads = 1;
    this->minQ = 0;
    this->profile = false;
    this->useQualities = false;
    this->splitReads = true;
    this->spanningReads = true;
    this->normalReads = true;
    this->mode = 2;
    this->gcCorrect = false;
    this->loadToMemory = false;
    this->sequenceDirectory = "";
    this->simCoverage = -1;
    this->stats = false;
}

ProgramOptions::ProgramOptions(
    std::string inFile, std::string fileList,
    std::string outFile, std::string vcfFile, std::string variantFile, std::string refFile, std::string sequenceDir, std::string samplingFile,
    bool wholeGenome, int nThreads, int minQ, int estimateDiffQ, int mode, bool verbose, bool outputDistributions, bool profile,
    bool noSplit, bool noSpanning, bool noStandard, bool noInsert, bool useQualities, bool gcCorrect, bool loadToMemory, int simCoverage, bool stats
)
{
    this->outFile = outFile;
    this->variantFile = variantFile;
    this->wholeGenome = wholeGenome;
    this->vcfFile = vcfFile;
    this->refFile = refFile;
    this->nThreads = std::max(nThreads, 1);
    this->verbose = verbose;
    this->outputDistributions = outputDistributions;
    this->profile = profile;
    this->estimateDiffQ = estimateDiffQ;
    this->noInsertSizes = noInsert;
    this->useQualities = useQualities;
    this->minQ = minQ;
    this->mode = mode;
    this->gcCorrect = gcCorrect;
    this->loadToMemory = loadToMemory;
    this->sequenceDirectory = sequenceDir;
    this->splitReads = !noSplit;
    this->spanningReads = !noSpanning;
    this->normalReads = !noStandard;
    this->simCoverage = simCoverage;
    this->stats = stats;

    determineBamFiles(inFile, fileList);
    determineSamplingRegions(samplingFile);
}

void ProgramOptions::determineBamFiles(std::string inFile, std::string fileList)
{
    if (inFile != "") {
        this->bamFileNames.push_back(inFile);
    } else {
        std::ifstream f(fileList);
        if (!f.is_open())
        {
            if (this->refFile != "")
                return;
            else
                throw std::runtime_error("Could not open list of BAM files provided.");
        }
        std::string filename;
        while (f.peek() != EOF)
        {
            std::getline(f, filename);
            this->bamFileNames.push_back(filename);
        }
        f.close();
    }
}

void ProgramOptions::determineSamplingRegions(std::string filename)
{
    if (filename == "")
        return;
    std::ifstream f(filename);
    if (!f.is_open())
        throw std::runtime_error("Could not open file with sampling regions.");
    std::string regionString;
    while (f.peek() != EOF)
    {
        std::getline(f, regionString);
        this->samplingRegions.push_back(GenomicRegion(regionString));
    }
    f.close();
}

std::vector<std::string> ProgramOptions::getBamFileNames()
{
    return this->bamFileNames;
}

std::string ProgramOptions::getOutputFileName()
{
    return this->outFile;
}

std::string ProgramOptions::getVariantFileName()
{
    return this->variantFile;
}

int ProgramOptions::getNumberOfThreads()
{
    return this->nThreads;
}

int ProgramOptions::getMinimumMappingQuality()
{
    return this->minQ;
}

bool ProgramOptions::isOptionVerbose()
{
    return this->verbose;
}

bool ProgramOptions::isOptionWriteUsedReads()
{
    return this->output;
}

bool ProgramOptions::isOptionUseWholeGenome()
{
    return this->wholeGenome;
}

bool ProgramOptions::isOptionWriteToFile()
{
    return this->outFile != "";
}

bool ProgramOptions::isOptionOutputDistributions()
{
    return this->outputDistributions;
}

bool ProgramOptions::isOptionNoSplitReads()
{
    return !this->splitReads;
}

bool ProgramOptions::isOptionNoSpanningReads()
{
    return !this->spanningReads;
}

bool ProgramOptions::isOptionNoNormalReads()
{
    return !this->normalReads;
}

std::string ProgramOptions::getVCFName()
{
    return this->vcfFile;
}

bool ProgramOptions::isOptionProfile()
{
    return this->profile;
}

std::string ProgramOptions::getRefFileName()
{
    return this->refFile;
}

bool ProgramOptions::isOptionGCCorrect()
{
    return this->gcCorrect;
}

int ProgramOptions::getEstimateDiffQual()
{
    return this->estimateDiffQ;
}

bool ProgramOptions::isOptionNoInsertSizes()
{
    return this->noInsertSizes;
}

bool ProgramOptions::isOptionUseQualities()
{
    return this->useQualities;
}

int ProgramOptions::getDistributionMode()
{
    return this->mode;
}

bool ProgramOptions::isOptionLoadToMemory()
{
    return this->loadToMemory;
}

void ProgramOptions::determineNumThreads(int nSamples, int nVariants)
{
    if (nVariants >= this->nThreads)
    {
        this->Tv = std::min(nVariants, this->nThreads);
        this->Ts = std::max(1, this->nThreads / this->Tv);
    } else {
        this->Ts = std::min(nSamples, this->nThreads);
        this->Tv = std::max(1, this->nThreads / this->Ts);
    }
}

int ProgramOptions::getNumberOfSampleThreads()
{
    return this->Ts;
}

int ProgramOptions::getNumberOfVariantThreads()
{
    return this->Tv;
}

std::string ProgramOptions::getSequenceDirectory()
{
    return this->sequenceDirectory;
}

std::vector<GenomicRegion> ProgramOptions::getSamplingRegions()
{
    return this->samplingRegions;
}

bool ProgramOptions::isOptionSample()
{
    return !(this->simCoverage < 0);
}

int ProgramOptions::getCoverage()
{
    return this->simCoverage;
}

bool ProgramOptions::isOptionStats()
{
    return this->stats;
}