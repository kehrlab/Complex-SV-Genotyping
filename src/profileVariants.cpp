#include "profileVariants.hpp"
#include "custom_types.hpp"
#include "profileHandler.hpp"
#include "seqan/arg_parse/argument_parser.h"
#include "variant.hpp"
#include <atomic>
#include <thread>
#include <unistd.h>

#ifndef DATE
#define DATE "1.1.1970"
#endif

#ifndef VERSION
#define VERSION "0.0"
#endif

int profileVariants(int argc, const char **argv)
{
    time_t now;
    std::string date;

    seqan::ArgumentParser argParser;
    seqan::ArgumentParser::ParseResult result = parseVariantProfileArgs(argParser, argc, argv);
    
    if (result == seqan::ArgumentParser::PARSE_HELP || 
        result == seqan::ArgumentParser::PARSE_VERSION ||
        result == seqan::ArgumentParser::PARSE_WRITE_CTD ||
        result == seqan::ArgumentParser::PARSE_EXPORT_HELP
        )
        return 0;
    else if (result != seqan::ArgumentParser::PARSE_OK)
        return 1;

    variantProfileParams params = getVariantProfileParameters(argParser);
    
    now = time(0);
    date = std::string(ctime(&now));
    date[date.find_last_of("\n")] = '\t';
    std::cout << date << "Gathering required parameters..." << std::endl;

    // get sMin, sMax and readLength from sample profiles
    int sMin {1000}, sMax {0}, readLength {-1};
    bool readLengthError {false};
    float insertMean = 0;
    std::vector<std::unordered_map<std::string, int32_t>> contigInfos;

    std::ifstream stream(params.sampleProfileFile);
    if (stream.is_open())
    {
        std::string filename;
        int counter = 0;
        while (stream.peek() != EOF)
        {
            ++counter;
            std::getline(stream, filename);

            Sample s;
            
            if (params.popdel) // initialize sample from popdel profile
            {
                PopDelProfileHandler popdelHandler(filename);
                s.loadPopDelProfile(popdelHandler);
                popdelHandler.close();
            } 
            else // read sample profile from profile-samples output
            {
                s.readSampleProfile(filename);
            }
            
            sMin = std::min(sMin, s.getLibraryDistribution().getMinInsert());
            sMax = std::max(sMax, s.getLibraryDistribution().getMaxInsert());
            insertMean += s.getLibraryDistribution().getInsertMean();
            contigInfos.push_back(s.getContigLengths());

            if (readLength < 0)
                readLength = s.getMaxReadLength();
            else if (readLength != s.getMaxReadLength())
            {
                readLengthError = true;
                std::cerr << "Read Length error in sample " << s.getSampleName() << ": " << s.getMaxReadLength() << std::endl;
            }
            s.close();
        }
        if (counter > 0)
            insertMean /= counter;
        stream.close();
    } else
    {
        std::cerr << "Could not open sample profile list " << params.sampleProfileFile << "." << std::endl;
        throw std::runtime_error("Aborted. Calculation of variant profiles requires access to sample profiles.");
    }

    if (readLengthError) {
	    std::cerr << "Consensus read length: " << readLength << std::endl;
	    throw std::runtime_error("Read lengths of sample profiles do not match. Cannot create common variant profiles.");
    }


    if (params.margin <= 0 && insertMean > 0)
        params.margin = insertMean;
    else if (params.margin <= 0)
        throw std::runtime_error("Error: Parameter -m missing and unable to determine m from sample profiles.");

    ContigInfo cInfo = mergeContigLengths(contigInfos);

    now = time(0);
    date = std::string(ctime(&now));
    date[date.find_last_of("\n")] = '\t';
    std::cout << date << "Loading variant descriptions and creating structures..." << std::endl;

    // load variant descriptions
    std::vector<complexVariant> variants;

    // for each given variant file
    for (std::string filename : params.variantFileNames)
    {
        // load variant descriptions from json file
        variantParser vParser(filename);
        std::vector<variantData> allVariantJunctions = vParser.getVariantJunctions();
        std::vector<std::string> variantNames = vParser.getVariantNames();
        std::vector<std::vector<std::string>> alleleNames = vParser.getAlleleNames();

        for (uint32_t i = 0; i < variantNames.size(); ++i)
        {
            variants.push_back(complexVariant(variantNames[i], alleleNames[i], allVariantJunctions[i], filename));
            variants[i].setFilterMargin(params.margin);
        }
    }

    // sort variants by number of novel junctions, making sure that larger variants are distributed as evenly as possible
    std::sort(variants.begin(), variants.end(), complexVariant::compareVariants); 


    // open output file for profile paths
    std::ofstream outStream(params.outFile);
    if (!outStream.is_open())
    {
        std::string msg = "Error: Could not open file " + params.outFile + " for writing.";
        throw std::runtime_error(msg.c_str());
    }

    // create variant profiles
    now = time(0);
    date = std::string(ctime(&now));
    date[date.find_last_of("\n")] = '\t';
    std::cout << date << "Calculating variant profiles..." << std::endl;

    params.nThreads = std::min((size_t) params.nThreads, variants.size()); // make sure not to create more threads than required

    std::vector<std::thread> threads(params.nThreads); // create pool of threads
    std::vector<bool> availability(params.nThreads, true); // track the availability of each thread
    std::atomic_int64_t counter = 0; // shared counter for tracking progress
    
    int overlap = 20;
    std::thread progressReporter(trackProgress, std::cref(counter), variants.size()); // one thread to track them all

    // Only handle one variant per thread instead of batch processing.
    // Memory allocated to each thread will be released once it is destroyed, 
    // avoiding build-up of unused memory that is not available to the remaining system.

    int64_t nextIdx = 0; // keep track of the next variant to be profiled
    for (int i = 0; i < params.nThreads; ++i) // first batch of threads is handled separately
    {
        availability[i] = false;
        threads[i] = std::thread(createVariantProfile, variants[i], params, overlap, readLength, sMin, sMax, cInfo, std::ref(counter), std::ref(availability), i);
        outStream << params.outDir + "/" + variants[i].getName() + ".profile" << std::endl;
        ++nextIdx;
    }

    while (nextIdx < variants.size())
    {
        for (int i = 0; i < params.nThreads; ++i) // cycle through all threads and assign the next variant to the next available thread 
        {
            if (availability[i]) {                                              // each thread only writes to i-th element -> no race conditions
                threads[i].join();                                              // if the thread is finished, we can join it immediately
                availability[i] = false;                                        // set availability BEFORE the new thread is started
                threads[i] = std::thread(createVariantProfile, variants[nextIdx], params, overlap, readLength, sMin, sMax, cInfo, std::ref(counter), std::ref(availability), i); 
                outStream << params.outDir + "/" + variants[nextIdx].getName() + ".profile" << std::endl; //  // write profile path to file 
                ++nextIdx;                                                      // next variant
                if (nextIdx == variants.size())
                    break;                                              
            }
        }

        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // check for available threads periodically
    }

    for (int i = 0; i < params.nThreads; ++i) // re-join all threads
        threads[i].join();
    progressReporter.join();


    outStream.close();


    now = time(0);
    date = std::string(ctime(&now));
    date[date.find_last_of("\n")] = '\t';
    std::cout << date << "Done" << std::endl;

    return 0;
}


seqan::ArgumentParser::ParseResult parseVariantProfileArgs(seqan::ArgumentParser &argParser, int argc, const char **argv)
{
    // arguments
    seqan::addArgument(argParser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "VARIANT FILE (list)")
    );
    seqan::addArgument(argParser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::OUTPUT_FILE, "OUTPUT FILE")
    );
    seqan::addArgument(argParser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::OUTPUT_DIRECTORY, "OUTPUT DIR")
    );
    seqan::addArgument(argParser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "SAMPLE_PROFILES")
    );

    seqan::setValidValues(argParser, 0, "JSON json txt");
    seqan::setValidValues(argParser, 1, "txt");
    seqan::setValidValues(argParser, 3, "txt");

    // options
    seqan::addOption(argParser, seqan::ArgParseOption(
        "m", "filter-margin", "Margin placed around reference breakpoints to obtain regions of interest. Default: 500.",
        seqan::ArgParseOption::ArgumentType::INTEGER, "MARGIN"
    ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "T", "threads", "Number of threads used for parallelization of variant profile calculation. Default: 1.",
        seqan::ArgParseOption::ArgumentType::INTEGER, "THREADS"
    ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "p", "popdel", "Use popdel profiles instead of sample profiles generated by 'ggtyper profile-samples'."
    ));

    // description
    seqan::addDescription(argParser, "Create profiles of given variants and write them to disk.");
    seqan::addUsageLine(argParser, "VARIANT_FILE OUTPUT_FILE OUTPUT_DIRECTORY SAMPLE_PROFILES [\033[4mOPTIONS\033[0m].");
    seqan::setDate(argParser, DATE);
    seqan::setVersion(argParser, VERSION);
    
    return seqan::parse(argParser, argc, argv);
}

variantProfileParams getVariantProfileParameters(const seqan::ArgumentParser &argParser)
{
    variantProfileParams params;

    std::string variantFileName;
    seqan::getArgumentValue(variantFileName, argParser, 0);

    if (seqan::getArgumentFileExtension(argParser, 0) == "txt")
    {
        std::ifstream stream(variantFileName);
        if (!stream.is_open())
        {
            std::cerr << "Could not open file " << variantFileName << " for reading." << std::endl;
            std::exit(1);
        }

        std::string filename;
        while (stream.peek() != EOF)
        {
            std::getline(stream, filename);
            params.variantFileNames.push_back(filename);
        }
        stream.close();
    } 
    else if (seqan::getArgumentFileExtension(argParser, 0) == "json")
    {
        params.variantFileNames.push_back(variantFileName);
    } 

    seqan::getArgumentValue(params.outFile, argParser, 1);
    seqan::getArgumentValue(params.outDir, argParser, 2);
    seqan::getArgumentValue(params.sampleProfileFile, argParser, 3);

    if (!std::filesystem::exists(params.outDir))
    {
        std::string msg = "ERROR: Directory '" + params.outDir + "' does not exist.";
        std::cerr << msg << std::endl;
        std::exit(1);
    }

    seqan::getOptionValue(params.margin, argParser, "filter-margin");
    seqan::getOptionValue(params.nThreads, argParser, "threads");
    seqan::getOptionValue(params.popdel, argParser, "popdel");
    return params;
}

ContigInfo mergeContigLengths(const std::vector<std::unordered_map<std::string, int32_t>> & contigInfos)
{
    std::unordered_set<std::string> cNames;
    for (auto & s : contigInfos)
        for (auto & chr : s)
            cNames.insert(chr.first);

    ContigInfo cInfo;

    for (auto & chr : cNames)
    {
        std::vector<int32_t> sizes;
        for (auto s : contigInfos)
            if (s.find(chr) != s.end())
                sizes.push_back(s[chr]);

        if (sizes.size() == contigInfos.size())
        {
            int32_t l = sizes[0];
            bool allMatch {true};
            for (uint32_t i = 1; i < sizes.size(); ++i)
            {
                if (sizes[i] != l)
                {
                    allMatch = false;
                    break;
                }
            }
            if (!allMatch)
                break;
            cInfo.cNames.push_back(chr);
            cInfo.cLengths.push_back(l);
        }
    }
    cInfo.sortNames();
    cInfo.calculateGlobalContigPositions();
    return cInfo;
}

inline void createVariantProfile(const complexVariant & variant, const variantProfileParams & params, int overlap, int readLength, int sMin, int sMax, ContigInfo cInfo, std::atomic_int64_t & counter, std::vector<bool> & availability, int idx)
{
    VariantProfile profile(variant, params.margin, overlap, readLength, sMin, sMax, cInfo);
    std::string profilePath = params.outDir + "/" + variant.getName() + ".profile";

    profile.calculateAlleleMasks();
    profile.writeProfile(profilePath);

    counter++;
    availability[idx] = true;
}

inline void trackProgress(const std::atomic_int64_t &counter, int64_t targetSize)
{
    time_t now;
    std::string date;

    while (counter < targetSize)
    {
        sleep(3);
        now = time(0);
        date = std::string(ctime(&now));
        date[date.find_last_of("\n")] = '\t';
        std::cout << "\r\e[K" << std::flush;
        std::cout << date << "Progress: " << counter << " / " << targetSize << " (";
        printf("%.2f", 100. * ((float) counter / targetSize));
        std::cout << "%)         " << std::flush;
    }
}