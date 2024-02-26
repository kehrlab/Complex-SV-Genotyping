#include "profileVariants.hpp"
#include "seqan/arg_parse/argument_parser.h"

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
    std::cout << date << "Loading variant descriptions and creating structures..." << std::endl;

    // load variants
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
    now = time(0);
    date = std::string(ctime(&now));
    date[date.find_last_of("\n")] = '\t';
    std::cout << date << "Gathering required parameters..." << std::endl;

    // get sMin, sMax and readLength from sample profiles
    int sMin {1000}, sMax {0}, readLength {-1};
    bool readLengthError {false};
    float insertMean = 0;
    std::vector<std::unordered_map<std::string, int>> contigInfos;

    if (params.sampleProfileFile != "")
    {
        std::ifstream stream(params.sampleProfileFile);
        if (stream.is_open())
        {
            std::string filename;
            int counter = 0;
            while (stream.peek() != EOF)
            {
                ++counter;
                std::getline(stream, filename);
                Sample s(filename);
                
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
        }
    } else {
        sMin = params.sMin;
        sMax = params.sMax;
        readLength = params.readLength;
    }
    if (readLengthError) {
	std::cerr << "Consensus read length: " << readLength << std::endl;
	throw std::runtime_error("Read lengths of sample profiles do not match. Cannot create common variant profiles.");
    }


    if (params.margin <= 0 && insertMean > 0)
        params.margin = insertMean;
    else if (params.margin <= 0)
        throw std::runtime_error("Error: Parameter -m missing and unable to determine m from sample profiles.");

    // create variant profiles
    now = time(0);
    date = std::string(ctime(&now));
    date[date.find_last_of("\n")] = '\t';
    std::cout << date << "Calculating variant profiles..." << std::endl;

    std::vector<std::string> profilePaths(variants.size());
    std::unordered_map<std::string, int> contigLengths = mergeContigLengths(contigInfos);
    int overlap = 20;

    #pragma omp parallel for num_threads(params.nThreads)
    for (uint32_t i = 0; i < variants.size(); ++i) 
    {
        VariantProfile profile(
            variants[i], params.margin, overlap, 
            readLength, sMin, sMax, contigLengths
            );
        profile.calculateAlleleMasks();
        profilePaths[i] = params.outDir + "/" + variants[i].getName() + ".profile";
        profile.writeProfile(profilePaths[i]);
    }

    // write profile paths to file
    std::ofstream stream(params.outFile);
    if (!stream.is_open())
    {
        std::string msg = "Error: Could not open file " + params.outFile + " for writing.";
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

    seqan::setValidValues(argParser, 0, "JSON json txt");
    seqan::setValidValues(argParser, 1, "txt");

    // options
    seqan::addOption(argParser, seqan::ArgParseOption(
        "S", "samples", 
        "List of sample profiles. Can be used to determine sMin, sMax and readLength. Otherwise these must be supplied manually.",
        seqan::ArgParseOption::ArgumentType::INPUT_FILE, "PROFILES"));
    seqan::setValidValues(argParser, "samples", "txt");

    seqan::addOption(argParser, seqan::ArgParseOption(
        "sMin", "min-insert-size", "Minimum insert size considered in profile creation. Default: 1.",
        seqan::ArgParseOption::ArgumentType::INTEGER, "SMIN"
    ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "sMax", "max-insert-size", "Maximum insert size considered in profile creation. Default: 1000.",
        seqan::ArgParseOption::ArgumentType::INTEGER, "SMAX"
    ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "l", "read-length", "Read length used to calculate profiles. Must match the samples that the profiles are genotyped later. Default: 100.",
        seqan::ArgParseOption::ArgumentType::INTEGER, "READLENGTH"
    ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "m", "filter-margin", "Margin placed around reference breakpoints to obtain regions of interest. Default: 500.",
        seqan::ArgParseOption::ArgumentType::INTEGER, "MARGIN"
    ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "T", "threads", "Number of threads used for parallelization of variant profile calculation. Default: 1.",
        seqan::ArgParseOption::ArgumentType::INTEGER, "THREADS"
    ));

    seqan::addDescription(argParser, "Create profiles of given variants and write them to disk.");
    seqan::addUsageLine(argParser, "VARIANT_FILE OUTPUT_FILE OUTPUT_DIRECTORY [\033[4mOPTIONS\033[0m].");
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

    if (!std::filesystem::exists(params.outDir))
    {
        std::string msg = "ERROR: Directory '" + params.outDir + "' does not exist.";
        std::cerr << msg << std::endl;
        std::exit(1);
    }

    if (! (seqan::isSet(argParser, "samples") || 
          (seqan::isSet(argParser, "min-insert-size") && 
           seqan::isSet(argParser, "max-insert-size") && 
           seqan::isSet(argParser, "read-length"))
          )
        )
    {
        throw std::runtime_error("Error in input parameters: Need to specify either '--samples' or '--min-insert-size', '--max-insert-size' and '--read-length'.");
    }

    seqan::getOptionValue(params.sampleProfileFile, argParser, "samples");
    seqan::getOptionValue(params.sMin, argParser, "min-insert-size");
    seqan::getOptionValue(params.sMax, argParser, "max-insert-size");
    seqan::getOptionValue(params.readLength, argParser, "read-length");
    seqan::getOptionValue(params.margin, argParser, "filter-margin");
    seqan::getOptionValue(params.nThreads, argParser, "threads");
    return params;
}

std::unordered_map<std::string, int> mergeContigLengths(
    const std::vector<std::unordered_map<std::string, int>> & contigInfos
    )
{
    std::unordered_set<std::string> cNames;
    for (auto & s : contigInfos)
        for (auto & chr : s)
            cNames.insert(chr.first);

    std::vector<std::string> cN;
    std::vector<int> cLengths;

    for (auto & chr : cNames)
    {
        std::vector<int> sizes;
        for (auto s : contigInfos)
        {
            if (s.find(chr) != s.end())
            {
                sizes.push_back(s[chr]);
            }
        }

        if (sizes.size() == contigInfos.size())
        {
            int l = sizes[0];
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
            cN.push_back(chr);
            cLengths.push_back(l);
        }
    }
    std::unordered_map<std::string, int> contigLengths;
    for (uint32_t i = 0; i < cN.size(); ++i)
        contigLengths[cN[i]] = cLengths[i];
    return contigLengths;
}
