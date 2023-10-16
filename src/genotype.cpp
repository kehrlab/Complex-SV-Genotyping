#include "genotype.hpp"
#include "genotypeResult.hpp"
#include "seqan/arg_parse/arg_parse_argument.h"
#include "seqan/arg_parse/arg_parse_option.h"
#include "variantProfile.hpp"
#include "vcfWriter.hpp"
#include <chrono>
#include <filesystem>
#include <stdexcept>

#ifndef DATE
#define DATE "1.1.1970"
#endif

#ifndef VERSION
#define VERSION "0.0"
#endif

int genotype(int argc, const char **argv)
{
    time_t now;
    std::string date;

    // create parser and get parameters
    seqan::ArgumentParser argParser;
    seqan::ArgumentParser::ParseResult result = parseGenotypeArgs(argParser, argc, argv);
    
    if (result == seqan::ArgumentParser::PARSE_HELP || 
        result == seqan::ArgumentParser::PARSE_VERSION ||
        result == seqan::ArgumentParser::PARSE_WRITE_CTD ||
        result == seqan::ArgumentParser::PARSE_EXPORT_HELP
        )
        return 0;
    else if (result != seqan::ArgumentParser::PARSE_OK)
        return 1;

    genotypeParameters params = getGenotypeParameters(argParser);


    now = time(0);
    date = std::string(ctime(&now));
    date[date.find_last_of("\n")] = '\t';
    std::cout << date << "Get sample parameters..." << std::endl;
    // load sample profile paths and check parameter consistency (sMin, sMax, readLength)
    std::vector<std::string> sampleProfiles;
    int sMin {-1}, sMax{-1}, readLength{-1};
    checkSampleParameters(sampleProfiles, sMin, sMax, readLength, params);

    // get list of variant profile files
    std::vector<std::string> variantFiles;
    std::string filename;
    std::ifstream vStream(params.variantList);
    if (!vStream.is_open())
        throw std::runtime_error("Could not open list of variant profiles for reading.");
    while (vStream.peek() != EOF)
    {
        std::getline(vStream, filename);
        variantFiles.push_back(filename);
    }
    vStream.close();
    
    // status
    now = time(0);
    date = std::string(ctime(&now));
    date[date.find_last_of("\n")] = '\t';
    std::cout << date << "Genotyping " << variantFiles.size() << " variants in ";
    std::cout << sampleProfiles.size() << " samples using " << params.nThreads << " threads..." << std::endl;

    // open output file and write header
    std::ofstream outFile(params.outputPrefix + "_genotype_results.tsv");
    if (!outFile.is_open())
    {
        std::string msg = "Could not open output file " + params.outputPrefix + "_genotype_results.tsv for writing.";
        throw std::runtime_error(msg.c_str());
    }
    outFile << "Variant\tSample\tFile\tGenotype\tMean_Quality\tLower_Bound\tUpper_Bound\tCertainty\tReads\tAvgMapQ\tMinMapQ\tMaxMapQ\tQualityPass";
    if (params.difficulties)
        outFile << "\tDifficulty";
    outFile << std::endl;

    // open VCF output
    VcfWriter vcfFile;
    if (params.vcfFile != "")
        vcfFile.setFileName(params.vcfFile);

    // genotyping
    int block {0};
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    while (block * params.nThreads < variantFiles.size())
    {
        // load variant profiles (currently all) and check parameters
        std::vector<VariantProfile> variantProfiles;
        loadVariantProfiles(variantProfiles, params, variantFiles, block);
        checkProfileParameters(sMin, sMax, readLength, variantProfiles, params);

        // temporary genotype results, required for VCF
        std::vector<std::vector<GenotypeResult>> tempResults;
        if (params.vcfFile != "")
            tempResults = std::vector<std::vector<GenotypeResult>>(
                variantProfiles.size(), 
                std::vector<GenotypeResult>(sampleProfiles.size())
                );

        // genotype current variants
        #pragma omp parallel for schedule(dynamic) collapse(2) num_threads(params.nThreads)
        for (int i = 0; i < variantProfiles.size(); ++i)
        {
            // in all samples
            for (int j = 0; j < sampleProfiles.size(); ++j)
            {
                Sample s(sampleProfiles[j]);
                GenotypeResult result(s.getFileName(), s.getSampleName(), s.getContigInfo().cNames, false);

                // load the relevant read pairs
                RecordManager bamRecords;
                BamFileHandler bamFileHandler(s.getFileName(), params.minMapQ);
                loadReadPairs(variantProfiles[i], bamRecords, bamFileHandler, s);

                // create the genotype distributions
                std::vector<std::string> genotypeNames;
                std::vector<GenotypeDistribution> genotypeDistributions;
                ReadPairFilter filter = variantProfiles[i].getFilter();
                std::vector<std::string> cNames = bamFileHandler.getContigInfo().cNames;
                LibraryDistribution & sampleDistribution = s.getLibraryDistribution();
                createGenotypeDistributions(
                    variantProfiles[i], 
                    genotypeNames, 
                    genotypeDistributions, 
                    0.001, 
                    cNames,
                    sampleDistribution
                );
        
                // calculate genotype difficulties based on KLD
                float difficulty = 0;
                if (params.difficulties)
                    calculateDifficulty(difficulty, genotypeNames, genotypeDistributions);

                // calculate the genotype likelihoods
                calculateGenotypeLikelihoods(
                    variantProfiles[i], 
                    genotypeNames, 
                    genotypeDistributions, 
                    bamRecords, 
                    filter, 
                    result
                );


                // call the genotype
                result.callGenotype();

                // write distributions if required
                if (params.distributions)
                    writeInsertSizeDistributions(variantProfiles[i].getVariant(), result, genotypeNames, genotypeDistributions);

                result.clearData();

                if (params.vcfFile != "")
                    tempResults[i][j] = result;

                // clean up
                bamFileHandler.closeInputFile();
                s.close();

                // write the genotype call to output file
                std::string outString = variantProfiles[i].getVariant().getName() + "\t" + result.getOutputString();
                if (params.difficulties)
                    outString += ("\t" + std::to_string(difficulty));
                #pragma omp critical
                outFile << outString << std::endl;
            }
        }

        // add VCF records
        if (params.vcfFile != "")
            for (int i = 0; i < variantProfiles.size(); ++i)
                vcfFile.addVariantRecords(tempResults[i], variantProfiles[i].getVariant());


        // Status and ETA
        std::chrono::steady_clock::time_point current = std::chrono::steady_clock::now();
        auto tAvg = std::chrono::duration_cast<std::chrono::seconds>(current-begin).count() / std::min((block + 1) * params.nThreads, (int) variantFiles.size());

        now = time(0);
        date = std::string(ctime(&now));
        date[date.find_last_of("\n")] = '\t';
        std::cout << "\r\e[K" << std::flush;
        std::cout << date << "Genotyped ";
        std::cout << std::min((block + 1) * params.nThreads, (int) variantFiles.size());
        std::cout << "/" << variantFiles.size() << " variants.\t\t\t";
        std::cout << "ETA: " << tAvg * (variantFiles.size() - std::min((block + 1) * params.nThreads, (int) variantFiles.size())) << "s     " << std::flush;

        // move to the next block of variants
        ++block;
    }
    std::cout << std::endl;
    outFile.close();

    // write VCF
    if (params.vcfFile != "")
    {
        vcfFile.write();
        vcfFile.closeVcfFile();
    }

    // write sample distributions
    if (params.distributions)
    {
        for (auto & p : sampleProfiles)
        {
            Sample s(p);
            s.getLibraryDistribution().writeDistribution(s.getFileName() + "_distributions/defaultDistribution.txt");
            s.close();
        }
    }
    
    return 0;
}

inline void loadVariantProfiles(std::vector<VariantProfile> & variantProfiles, genotypeParameters & params, std::vector<std::string> & profilePaths, int block)
{   
    int beginIdx = block * params.nThreads;
    int endIdx = std::min(((block + 1) * params.nThreads - 1), (int) profilePaths.size() - 1);
    std::vector<VariantProfile>().swap(variantProfiles);

    for (int i = beginIdx; i <= endIdx; ++i)
        variantProfiles.push_back(VariantProfile(profilePaths[i]));
    if (variantProfiles.size() == 0)
        throw std::runtime_error("No variant profiles successfully read.");
    return;
}

inline void checkProfileParameters(int & sMin, int & sMax, int & readLength, std::vector<VariantProfile> & variantProfiles, genotypeParameters & params)
{
    for (auto it = variantProfiles.begin(); it != variantProfiles.end(); )
    {
        if (it->getMinInsert() > sMin || it->getMaxInsert() < sMax)
        {
            std::string msg = "At least one sample contains insert sizes not covered by profile for variant " + it->getName() + ". Dropping variant.";
            std::cerr << msg << std::endl;
            it = variantProfiles.erase(it);
            continue;
        }
        if (it->getReadLength() != readLength)
        {
            std::string msg = "Read length of profile " + 
                it->getName() + 
                " does not match the given samples. Dropping variant. " +
                "Make sure that variant profiles are generated with parameters matching the profiles of samples on which they will be used.";
            std::cerr << msg << std::endl;
            it = variantProfiles.erase(it);
            continue;
        }
        if (!it->variantStructureIsPresent())
        {
            std::cout << "There is no structure present?" << std::endl;
            std::cout << it->getName() << "\t" << it->variantStructureIsPresent() << std::endl;
            if (params.variantFile == "")
            {
                std::string msg = "Variant description could not be loaded from profile. A variant description file needs to be specified.";
                throw std::runtime_error(msg.c_str());
            } else {
                if (!it->loadVariantStructure(params.variantFile, it->getName()))
                {
                    std::string msg =  "Description of variant " + it->getName() + " could not be loaded from file " + params.variantFile + ".";
                    throw std::runtime_error(msg.c_str());
                } else {
                    it->getFilter() = ReadPairFilter(it->getVariant().getAllBreakpoints(), it->getMargin(), 0);
                }
            }
        }
        ++it;
    }
    return;
}

inline void checkSampleParameters(std::vector<std::string> & sampleProfiles, int & sMin, int & sMax, int & readLength, genotypeParameters & params)
{
    std::string filename;
    std::ifstream stream(params.sampleList);
    if (!stream.is_open())
        throw std::runtime_error("Could not open list of sample profiles for reading.");
    while (stream.peek() != EOF)
    {
        std::getline(stream, filename);
        sampleProfiles.push_back(filename);
        Sample s(filename);
        if (readLength < 0)
            readLength = s.getLibraryDistribution().getReadLength();
        else if (s.getLibraryDistribution().getReadLength() != readLength)
            throw std::runtime_error("Read length in samples does not match. Variant profiles for these samples must be generated separately.");
        sMax = std::max(sMax, s.getLibraryDistribution().getMaxInsert());
        if (sMin < 0)
            sMin = s.getLibraryDistribution().getMinInsert();
        else
            sMin = std::min(s.getLibraryDistribution().getMinInsert(), sMin);
        s.close();
    }
    stream.close();
    return;
}

inline void loadReadPairs(VariantProfile & variantProfile, RecordManager & recordManager, BamFileHandler & bamFileHandler, Sample & s)
{
    bamFileHandler.setRegions(
        variantProfile.getVariant().calculateAssociatedRegions(s.getLibraryDistribution())
    );
    recordManager.setReadPairs(bamFileHandler.get_read_pairs());
    return;
}

inline void createGenotypeDistributions(VariantProfile & variantProfile, std::vector<std::string> & genotypeNames, std::vector<GenotypeDistribution> & genotypeDistributions, float eps, std::vector<std::string> & cNames, LibraryDistribution & sampleDistribution)
{
    std::unordered_map<std::string, GenotypeDistribution> gtDists;
    variantProfile.calculateGenotypeDistributions(gtDists, sampleDistribution, 0.001);
    
    for (auto & dist : gtDists)
    {
        dist.second.setPossibleContigs(cNames);
        genotypeNames.push_back(dist.first);
        genotypeDistributions.push_back(dist.second);
    }
    std::unordered_map<std::string, GenotypeDistribution>().swap(gtDists);
}

inline void calculateGenotypeLikelihoods(
    VariantProfile & variantProfile, std::vector<std::string> & genotypeNames, 
    std::vector<GenotypeDistribution> & genotypeDistributions, 
    RecordManager & bamRecords, ReadPairFilter & filter, GenotypeResult & result
    )
{
    result.initLikelihoods(genotypeNames);
    for (auto it = bamRecords.getTemplates().begin(); it != bamRecords.getTemplates().end(); ++it)
    {
        if (it->isInterestingReadPair(filter)) 
        {
            it->findSpanningReads(variantProfile.getVariant().getAllBreakpoints());
            it->findSplitReads(variantProfile.getVariant().getAllJunctions(), variantProfile.getChromosomeStructures());
            it->determineLocationStrings();
            adjustLikelihoods(*it, genotypeNames, genotypeDistributions, result);
        }
    }
}

inline void adjustLikelihoods(ReadTemplate & readTemplate, std::vector<std::string> & genotypeNames, std::vector<GenotypeDistribution> & genotypeDistributions, GenotypeResult & result)
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

    std::vector<float> probabilities;
    bool insertSizeWithinLimits = false;
    for (auto dist : genotypeDistributions) {
	if (insertSize >= dist.getMinInsertSize() - 500 && insertSize <= dist.getMaxInsertSize() + 500)
		insertSizeWithinLimits = true;
        probabilities.push_back(
            std::max(dist.getProbability(insertSize, orientation, split, spanning, interChromosome, regionString, junctionString, breakpointString, rNamePairs, true), minValue)
            );
    } 
    
    if (insertSizeWithinLimits)
    	result.storeEvidence(insertSize, orientation, split, spanning, interChromosome, regionString, junctionString, breakpointString, rNamePairs, mappingQualities);	
    result.addTemplateProbabilities(genotypeNames, probabilities, readTemplate.getTemplateWeight());

    return;
}

void calculateDifficulty(float & difficulty, std::vector<std::string> & genotypeNames, std::vector<GenotypeDistribution> & genotypeDistributions)
{
    float minKLD = 1000;
    for (int i = 0; i < genotypeDistributions.size(); ++i)
        for (int j = 0; j < genotypeDistributions.size(); ++j)
            if (j != i)
                minKLD = std::min(minKLD, genotypeDistributions[i].calculateKLD(genotypeDistributions[j]));
    difficulty = 1000 / minKLD;
}


void writeInsertSizeDistributions(complexVariant & variant, GenotypeResult & result, std::vector<std::string> & genotypeNames, std::vector<GenotypeDistribution> & genotypeDistributions)
{
    std::string variantName = "unnamed";
    if (variant.getName() != "")
        variantName = variant.getName();
    
    std::string directory = result.getFilename() + "_distributions/";
    try {
        std::filesystem::create_directory(directory);
    } catch (std::filesystem::filesystem_error)
    {
        std::cerr << "Could not create variant distribution directory" << std::endl;
        return;
    }

    directory += (variantName + "/");
    try {
        std::filesystem::create_directory(directory);
    }catch (std::filesystem::filesystem_error)
    {
        std::cerr << "Could not create variant distribution directory" << std::endl;
        return;
    }

    // write theoretical distributions
    for (int i = 0; i < genotypeNames.size(); ++i)
    {
        std::string prefix = directory + "numericalDistribution";
        
        std::string gtName = genotypeNames[i];
        int slashIndex = gtName.find("/");
        gtName.replace(slashIndex, 1, "-");

        if (variantName != "")
            prefix.append("_" + variantName);
        prefix.append("_" + gtName);
        genotypeDistributions[i].writeDistributionBinned(prefix);
    }

    // write observed distribution
    std::string observedDistributionPrefix = directory + "observedDistribution_" + variantName;
    result.writeEvidence(observedDistributionPrefix);

    // write bootstrap result
    std::string bootstrapPrefix = directory + variantName;
    result.writeBootstrapData(bootstrapPrefix);
}


seqan::ArgumentParser::ParseResult parseGenotypeArgs(seqan::ArgumentParser &argParser, int argc, const char **argv)
{
    seqan::addArgument(argParser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "variant profiles"));
    seqan::addArgument(argParser, seqan::ArgParseArgument(seqan::ArgParseArgument::INPUT_FILE, "sample profiles"));
    seqan::addArgument(argParser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "output prefix"));
    seqan::setValidValues(argParser, 0, "txt");
    seqan::setValidValues(argParser, 1, "txt");

    seqan::addOption(argParser, seqan::ArgParseOption(
        "T", "threads", 
        "Number of threads to use for paralleliation.",
        seqan::ArgParseOption::ArgumentType::INTEGER, "THREADS"
        ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "d", "write-distributions", 
        "Write theoretical and observed distributions to disk. Location of the input bam file must be writeable."
        ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "e", "estimate-difficulty", 
        "Give a difficulty estimation of the genotyped variants."
        ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "q", "minimum-quality", 
        "Required quality for alignments to be considered during genotyping. Default: 0.",
        seqan::ArgParseOption::ArgumentType::INTEGER, "MINMAPQ"
        ));
    seqan::addOption(argParser, seqan::ArgParseOption(
        "V", "vcf-out", "Path of desired vcf output file. Note: VCF output is experimental and incomplete, not recommended.",
        seqan::ArgParseOption::OUTPUT_FILE, "VCF"
    ));
    seqan::setValidValues(argParser, "vcf-out", "vcf VCF");

    seqan::addDescription(argParser, "Genotype given variants (specified by profiles) in all samples (specified by profiles).");
    seqan::addUsageLine(argParser, "ggtyper genotype VARIANT_PROFILES SAMPLE_PROFILES OUTPUT_PREFIX [\033[4mOPTIONS\033[0m].");
    seqan::setDate(argParser, DATE);
    seqan::setVersion(argParser, VERSION);
    
    return seqan::parse(argParser, argc, argv);
}


genotypeParameters getGenotypeParameters(seqan::ArgumentParser &argParser)
{
    genotypeParameters params;
    seqan::getArgumentValue(params.variantList, argParser, 0);
    seqan::getArgumentValue(params.sampleList, argParser, 1);
    seqan::getArgumentValue(params.outputPrefix, argParser, 2);

    seqan::getOptionValue(params.nThreads, argParser, "threads");
    seqan::getOptionValue(params.difficulties, argParser, "estimate-difficulty");
    seqan::getOptionValue(params.minMapQ, argParser, "minimum-quality");
    seqan::getOptionValue(params.distributions, argParser, "write-distributions");
    seqan::getOptionValue(params.vcfFile, argParser, "vcf-out");

    return params;
}
