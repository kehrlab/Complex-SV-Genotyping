#include "profileSamples.hpp"
#include "seqan/arg_parse/arg_parse_option.h"
#include "seqan/arg_parse/argument_parser.h"
#include <cstdlib>
#include <filesystem>
#include <stdexcept>

#ifndef DATE
#define DATE "1.1.1970"
#endif

#ifndef VERSION
#define VERSION "0.0"
#endif

int profileSamples(int argc, const char ** argv)
{
    time_t now;
    std::string date;

    // parser command line arguments
    seqan::ArgumentParser argParser;
    
    seqan::ArgumentParser::ParseResult result = parseSampleProfileArgs(argParser, argc, argv);

    if (result == seqan::ArgumentParser::PARSE_HELP || 
        result == seqan::ArgumentParser::PARSE_VERSION ||
        result == seqan::ArgumentParser::PARSE_WRITE_CTD ||
        result == seqan::ArgumentParser::PARSE_EXPORT_HELP
        )
        return 0;
    else if (result != seqan::ArgumentParser::PARSE_OK)
        return 1;
    
    sampleProfileParams params = getSampleProfileParameters(argParser);

    now = time(0);
    date = std::string(ctime(&now));
    date[date.find_last_of("\n")] = '\t';
    std::cout << date << "Gathering required information about BAM files and regions..." << std::endl;

    // get contig infos from bam files
    std::vector<ContigInfo> contigInfos;
    std::vector<int> readLengths(params.sampleFiles.size());
    std::vector<int> readGroups(params.sampleFiles.size());
    bool rgError {false};
    for (int i = 0; i < params.sampleFiles.size(); ++i)
    {
        BamFileHandler tempFile(params.sampleFiles[i]);
        readGroups[i] = tempFile.getRGNumber();
        readLengths[i] = tempFile.getReadLength();

        if (readGroups[i] != 1 && !params.forceRG)
        {
            rgError = true;
            std::cerr << "Number of read groups (" << readGroups[i];
            std::cerr << ") in sample " << params.sampleFiles[i] << " is not 1." << std::endl;
        }

        contigInfos.push_back(tempFile.getContigInfo());
        tempFile.closeInputFile();
    }

    if (rgError)
    {
        std::cerr << "Abort. Samples with more than one read group detected. This may lead to errors. To ignore this, use flag '-f'." << std::endl;
        return 1;
    }
    
    // determine regions from which insert size distribution is calclated
    std::vector<GenomicRegion> regions;
    RegionSampler sampler(contigInfos, params.regions);
    if (!params.wholeGenome)
        regions = sampler.sampledInsertRegions();
    else
        regions = params.regions;

    // calculate distribution for each file and write it to disk
    std::vector<std::string> profilePaths(params.sampleFiles.size());
    
    now = time(0);
    date = std::string(ctime(&now));
    date[date.find_last_of("\n")] = '\t';
    std::cout << date << "Creating sample profiles..." << std::endl;

    #pragma omp parallel for num_threads(params.nThreads)
    for (int i = 0; i < params.sampleFiles.size(); ++i)
    {
        Sample s(params.sampleFiles[i], regions);
        profilePaths[i] = params.outDir + "/" + s.getSampleName() + ".profile";
        s.writeSampleProfile(profilePaths[i]);
        s.close();
    }

    // write profile paths to disk
    std::ofstream stream(params.outFile);
    if (!stream.is_open())
    {
        std::string msg = "Could not open file " + params.outFile + " for writing.";
        throw std::runtime_error(msg.c_str());
    }
    for (auto & f : profilePaths)
        stream << f << std::endl;
    stream.close();

    now = time(0);
    date = std::string(ctime(&now));
    date[date.find_last_of("\n")] = '\t';
    std::cout << date << "Done" << std::endl;

    return 0;
}



seqan::ArgumentParser::ParseResult parseSampleProfileArgs(seqan::ArgumentParser &argParser, int argc, const char **argv)
{
    // add the required arguments
    seqan::addArgument(argParser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "bam file(list)"));
    seqan::addArgument(argParser, seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUT_FILE, "output file"));
    seqan::addArgument(argParser, seqan::ArgParseArgument(seqan::ArgParseArgument::OUTPUT_DIRECTORY, "output dir"));

    seqan::setValidValues(argParser, 0, "txt bam");
    seqan::setValidValues(argParser, 1, "txt");

    // add the options
    seqan::addOption(argParser, seqan::ArgParseOption(
        "w", "whole-genome", 
        "If specified, all records within the given bam file(s) or the regions specified by option --regions are used."
        ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "r", "regions", 
        "The regions to use for creation of the library distribution. Must be a '.txt' file with one region per line, e.g., chr1:10000-20000.",
        seqan::ArgParseOption::ArgumentType::INPUT_FILE, "REGIONS"
        ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "T", "threads", 
        "Number of threads used to parallelize sample profile calculation. Default: 1.",
        seqan::ArgParseOption::ArgumentType::INTEGER, "THREADS"
        ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "f", "force-rg",
        "Ignore presence of multiple read groups and force calculation of library distribution across all read groups."
    ));

    seqan::addDescription(argParser, "Create profiles of given bam files and write them to disk.");
    seqan::addUsageLine(argParser, "BAM_FILE(LIST) OUTPUT_FILE OUTPUT_DIR [\033[4mOPTIONS\033[0m].");
    seqan::setDate(argParser, DATE);
    seqan::setVersion(argParser, VERSION);

    return seqan::parse(argParser, argc, argv);;
}



sampleProfileParams getSampleProfileParameters(const seqan::ArgumentParser &argParser)
{
    sampleProfileParams params;

    // get list of bam files
    std::string bamSource;
    seqan::getArgumentValue(bamSource, argParser, 0);
    if (seqan::getArgumentFileExtension(argParser, 0) == "bam")
    {
        params.sampleFiles.push_back(bamSource);
    }
    else
    {
        std::string filename;
        std::ifstream stream(bamSource);
        if (!stream.is_open())
            throw std::runtime_error("Error: Could not open list of bam files for reading.");
        while (stream.peek() != EOF)
        {
            std::getline(stream, filename);
            params.sampleFiles.push_back(filename);
        }
        stream.close();
    }

    // get the output arguments
    seqan::getArgumentValue(params.outFile, argParser, 1);
    seqan::getArgumentValue(params.outDir, argParser, 2);
    if (!std::filesystem::exists(params.outDir))
    {
        std::string msg = "ERROR: Directory '" + params.outDir + "' does not exist.";
        std::cerr << msg << std::endl;
        std::exit(1);
    }
    
    // get whole genome flag
    seqan::getOptionValue(params.wholeGenome, argParser, "whole-genome");
    seqan::getOptionValue(params.nThreads, argParser, "threads");
    seqan::getOptionValue(params.forceRG, argParser, "force-rg");

    // get regions (if any)
    if (seqan::isSet(argParser, "regions"))
    {
        std::string regionFile;
        std::string regionString;
        seqan::getOptionValue(regionFile, argParser, "regions");
        std::ifstream stream(regionFile);
        if (!stream.is_open())
            throw std::runtime_error("Error: Could not open region file.");
        while (stream.peek() != EOF)
        {
            std::getline(stream, regionString);
            params.regions.push_back(GenomicRegion(regionString));
        }
        stream.close();
    }
    return params;
}
