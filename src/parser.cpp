#include "custom_types.hpp"

#include "options.hpp"
#include "parser.hpp"
#include "seqan/arg_parse/arg_parse_argument.h"
#include "seqan/arg_parse/arg_parse_option.h"
#include "seqan/arg_parse/argument_parser.h"
#include <stdexcept>

using namespace seqan;

parser::parser()
{

}

parser::parser(int argc, const char **argv)
{
    addOptionsIO();
    addOptionsParams();
    parseOptions(argc, argv);
    if (this->parsingStatus == ArgumentParser::PARSE_OK)
    {
        checkForArgumentConflict();
        extractOptions();
    }   
    checkParsingStatus();
    checkGCRequirements();
}

void 
parser::addOptionsIO()
{
    addOptionInputFile();
    addOptionInputFileList();
    addOptionOutputFile();
    addOptionConfigFile();
    addOptionOutputVCF();
    addOptionGenomeFile();
    addOptionSequenceDirectory();
    addOptionSamplingRegions();
}

void 
parser::addOptionsParams()
{
    addOptionNumThreads();
    addOptionVerbose();
    addOptionWholeGenome();
    addOptionOutputInsertSizeDistributions();
    addOptionProfiling();
    addOptionEstimateDifficulty();
    addOptionNoSplits();
    addOptionNoSpanning();
    addOptionNoStandard();
    addOptionNoInsertSize();
    //addOptionUseQualities();
    addOptionMinQuality();
    addOptionMode();
    addOptionGCCorrect();
    addOptionLoadToMemory();
    addOptionCoverage();
    addOptionStats();
}

void 
parser::addOptionInputFile()
{
    addOption(this->argParser, ArgParseOption(
        "i", "input-file", "BAM input file",
        ArgParseArgument::INPUT_FILE, "IN"
    ));
    setValidValues(this->argParser, "input-file", "bam");
}

void 
parser::addOptionInputFileList()
{
    addOption(this->argParser, ArgParseOption(
        "F", "file-list", "text fiel containing list of bam files",
        ArgParseOption::INPUT_FILE, "TXT"
    ));
    setValidValues(this->argParser, "file-list", "txt");
}

void 
parser::addOptionOutputFile()
{
    addOption(this->argParser, ArgParseOption(
        "o", "output-file", "output file for genotype information",
        ArgParseOption::OUTPUT_FILE, "OUTFILE"
    ));
}

void
parser::addOptionOutputVCF()
{
    addOption(this->argParser, ArgParseOption(
        "V", "output-vcf", "VCF output for variant breakends",
        ArgParseOption::OUTPUT_FILE, "VCF"
    ));
    setValidValues(this->argParser, "output-vcf", "vcf");
}

void parser::addOptionGenomeFile()
{
    addOption(this->argParser, ArgParseOption(
        "G", "reference-genome", "fasta file containing reference sequences",
        ArgParseOption::INPUT_FILE, "FASTA"
    ));
    setRequired(this->argParser, "reference-genome", true);
}


void 
parser::addOptionConfigFile()
{
    addOption(this->argParser, ArgParseOption(
        "C", "config-file", "json file containing novel junctions",
        ArgParseOption::INPUT_FILE, "JSON"
    ));
    setValidValues(this->argParser, "config-file", "json JSON");
}

void parser::addOptionSequenceDirectory()
{
    addOption(this->argParser, ArgParseOption(
        "s", "sequence-directory", "directory that created allele sequences are written to",
        ArgParseOption::OUTPUT_DIRECTORY, "DIR"
    ));
}

void parser::addOptionSamplingRegions()
{
    addOption(this->argParser, ArgParseOption(
        "R", "sampling-regions", "text file containing all regions used for sampling default insert sizes, one region per line",
        ArgParseOption::INPUT_FILE, "REGIONS"
    ));
}

void 
parser::addOptionNumThreads()
{
    addOption(this->argParser, ArgParseOption(
        "T", "num-threads", "max number of threads to use",
        ArgParseOption::INTEGER, "THREADS"
    ));
}

void parser::addOptionCoverage()
{
    addOption(this->argParser, ArgParseOption(
        "z", "coverage", "average sampling coverage for determination of theoretical distribution; if not given no sampling is performed",
        ArgParseOption::INTEGER, "COVERAGE"
    ));
}

void
parser::addOptionMinQuality()
{
    addOption(this->argParser, ArgParseOption(
        "minQ", "minMapQ", "minimum mapping quality used to filter records (inclusive)",
        ArgParseOption::INTEGER, "MAPQ" 
    ));
}

void parser::addOptionMode()
{
    addOption(this->argParser, ArgParseOption(
        "mode", "distribution-mode", "internal storing of insert size distributions (0: one distribution for split and spanning, each; 1: split and spanning parted by orientation; 2: like 0, but location is considered; 3: insert size is only used for RF, otherwise like 0)",
        ArgParseOption::INTEGER, "MODE"
    ));
}

void 
parser::addOptionVerbose()
{
    addOption(this->argParser, ArgParseOption(
        "v", "verbose", "print additional information"
    ));
}

void
parser::addOptionLoadToMemory()
{
    addOption(this->argParser, ArgParseOption(
        "M", "load-to-memory", "load reference sequences of sampled chromosomes into memory"
    ));
}

void parser::addOptionProfiling()
{
    addOption(this->argParser, ArgParseOption(
        "p", "profile", "profile performance"
    ));
}

void 
parser::addOptionWholeGenome()
{
    addOption(this->argParser, ArgParseOption(
        "w", "whole-genome", 
        "use all records instead of sampling to determine insert size distribution"
    ));
}

void parser::addOptionGCCorrect()
{
    addOption(this->argParser, ArgParseOption(
        "gc", "gc-correct", "perform gc bias correction based on observed (sampled) records ind the bam file"
    ));
}

void parser::addOptionOutputInsertSizeDistributions()
{
    addOption(this->argParser, ArgParseOption(
        "d", "output-distributions", "write insert size distributions for all genotypes to text files"
    ));
}

void parser::addOptionEstimateDifficulty()
{
    addOption(this->argParser, ArgParseOption(
        "e", "estimate-difficulty", "get an estimate of how many reads are needed to genotype variant with a certain quality",
        ArgParseOption::INTEGER, "DESIRED QUALITY"
    ));
}

void parser::addOptionNoSplits()
{
    addOption(this->argParser, ArgParseOption(
        "k", "no-split", "do not use split reads"
    ));
}

void parser::addOptionNoSpanning()
{
    addOption(this->argParser, ArgParseOption(
        "l", "no-spanning", "do not use split reads"
    ));
}

void parser::addOptionNoStandard()
{
    addOption(this->argParser, ArgParseOption(
        "m", "no-standard", "do not use split reads"
    ));
}

void parser::addOptionNoInsertSize()
{
    addOption(this->argParser, ArgParseOption(
        "n", "no-insert-sizes", "combine each insert size distribution into one probability for the corresponding read pair class"
    ));
}

void parser::addOptionUseQualities()
{
    addOption(this->argParser, ArgParseOption(
        "q", "use-qualities", "use alignment qualities during bootstrapping"
    ));
}

void parser::addOptionStats()
{
    addOption(this->argParser, ArgParseOption(
        "x", "stats", "record stats per read group"
    ));
}

void 
parser::parseOptions(int argc, const char **argv)
{
    this->parsingStatus = seqan::parse(this->argParser, argc, argv);
}

void parser::checkForArgumentConflict()
{
    if (singleFileAndListPresent() || singleFileAndListMissing())
    {
        this->parsingStatus = ArgumentParser::PARSE_ERROR;
    }
}

void parser::checkGCRequirements()
{
    if (isSet(this->argParser, "gc-correct") && !isSet(this->argParser, "reference-genome"))
    {
        std::cerr << "GC correction requires specification of reference genome (Option -G)" << std::endl;
        this->parsingStatus = ArgumentParser::PARSE_ERROR;
    }
}

bool parser::singleFileAndListPresent()
{
    if (isSet(this->argParser, "file-list") && isSet(this->argParser, "input-file"))
    {
        std::cerr << "Conflicting Arguments. Set either -i or -F, not both." << std::endl;
        return(true);
    }
    return false;
}

bool parser::singleFileAndListMissing()
{
    if (!(isSet(this->argParser, "file-list") || isSet(this->argParser, "input-file")))
    {
        if (! (isSet(this->argParser, "reference-genome"))) {
            std::cerr << "Missing Argument. Please set either -i or -F." << std::endl;
            return true;
        }
    }
    return false;
}

void parser::checkParsingStatus ()
{
    if (this->parsingStatus == ArgumentParser::PARSE_ERROR)
        throw std::runtime_error("Argument Parsing Failed. Please check for conflicting arguments and consult the help page.");
}

bool parser::wasSuccessful()
{
    return (this->parsingStatus == ArgumentParser::PARSE_OK);
}

void parser::extractOptions()
{
    std::string inFile, outFile, insertSizeFile, fileList, variantFile, vcfFile, refFile, sequenceDir, samplingFile;
    int nThreads = 1, minQ = 0, estimate = 0, mode = 2;
    bool wholeGenome = false; 
    bool verbose = false, outputDistributions = false, profile = false;
    bool noSplit = false, noSpanning = false, noStandard = false, noInsert = false;
    bool useQualities = false;
    bool gcCorrect = false;
    bool loadToMemory = false;
    int coverage = -1;
    bool stats = false;

    getOptionValue(inFile, this->argParser, "input-file");
    getOptionValue(outFile, this->argParser, "output-file");
    getOptionValue(variantFile, this->argParser, "config-file");
    getOptionValue(fileList, this->argParser, "file-list");
    getOptionValue(nThreads, this->argParser, "num-threads");
    getOptionValue(wholeGenome, this->argParser, "whole-genome");
    getOptionValue(verbose, this->argParser, "verbose");
    getOptionValue(outputDistributions, this->argParser, "output-distributions");
    getOptionValue(vcfFile, this->argParser, "output-vcf");
    getOptionValue(profile, this->argParser, "profile");
    getOptionValue(refFile, this->argParser, "reference-genome");
    getOptionValue(gcCorrect, this->argParser, "gc-correct");
    getOptionValue(estimate, this->argParser, "estimate-difficulty");
    getOptionValue(noSplit, this->argParser, "no-split");
    getOptionValue(noSpanning, this->argParser, "no-spanning");
    getOptionValue(noStandard, this->argParser, "no-standard");
    getOptionValue(noInsert, this->argParser, "no-insert-sizes");
    //getOptionValue(useQualities, this->argParser, "use-qualities");
    getOptionValue(minQ, this->argParser, "minQ");
    getOptionValue(mode, this->argParser, "distribution-mode");
    getOptionValue(loadToMemory, this->argParser, "load-to-memory");
    getOptionValue(sequenceDir, this->argParser, "sequence-directory");
    getOptionValue(samplingFile, this->argParser, "sampling-regions");
    getOptionValue(coverage, this->argParser, "coverage");
    getOptionValue(stats, this->argParser, "stats");
    
    this->options = ProgramOptions(
        inFile, fileList, outFile, vcfFile,
        variantFile, refFile,
        sequenceDir, samplingFile,
        wholeGenome, nThreads, minQ, estimate, mode, verbose,
        outputDistributions, profile,
        noSplit, noSpanning, noStandard, noInsert, 
        useQualities, gcCorrect, loadToMemory, coverage, stats
    );
}

ProgramOptions parser::getOptions()
{
    return this->options;
}