#include "genotype.hpp"

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
    std::cout << date << "\tLoading variant profiles..." << std::endl;

    // load variant profiles (currently all)
    std::vector<VariantProfile> variantProfiles;
    loadVariantProfiles(variantProfiles, params);
    
    // check consistency of profile parameters (readLength, sMin, sMax)
    int sMinVar {0}, sMaxVar {0}, readLengthVar {0};
    checkProfileParameters(sMinVar, sMaxVar, readLengthVar, variantProfiles, params);

    now = time(0);
    date = std::string(ctime(&now));
    std::cout << date << "\tGet sample parameters..." << std::endl;

    // load sample profile paths and check parameter consistency (sMin, sMax, readLength)
    std::vector<std::string> sampleProfiles;
    checkSampleParameters(sampleProfiles, sMinVar, sMaxVar, readLengthVar, params);
    

    now = time(0);
    date = std::string(ctime(&now));
    std::cout << date << "\tGenotyping " << variantProfiles.size() << " variants in ";
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

    // genotype all variants
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
            {
                writeInsertSizeDistributions(variantProfiles[i].getVariant(), result, genotypeNames, genotypeDistributions);
                sampleDistribution.writeDistribution(result.getFilename() + "_distributions/defaultDistribution.txt");
            }

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
    outFile.close();

    now = time(0);
    date = std::string(ctime(&now));
    std::cout << date << "\tDone" << std::endl;

    return 0;
}

inline void loadVariantProfiles(std::vector<VariantProfile> & variantProfiles, genotypeParameters & params)
{
    std::string filename;
    std::ifstream vStream(params.variantList);
    if (!vStream.is_open())
        throw std::runtime_error("Could not open list of variant profiles for reading.");
    while (vStream.peek() != EOF)
    {
        std::getline(vStream, filename);
        variantProfiles.push_back(VariantProfile(filename));
    }
    vStream.close();
    if (variantProfiles.size() == 0)
        throw std::runtime_error("No variant profiles successfully read.");
    return;
}

inline void checkProfileParameters(int & sMinVar, int & sMaxVar, int & readLengthVar, std::vector<VariantProfile> & variantProfiles, genotypeParameters & params)
{
    sMinVar = variantProfiles[0].getMinInsert();
    sMaxVar = variantProfiles[0].getMaxInsert();
    readLengthVar = variantProfiles[0].getReadLength();

    for (int i = 1; i < variantProfiles.size(); ++i)
    {
        sMinVar = std::max(sMinVar, variantProfiles[i].getMinInsert());
        sMaxVar = std::min(sMaxVar, variantProfiles[i].getMaxInsert());
        if (variantProfiles[i].getReadLength() != readLengthVar)
            throw std::runtime_error(
                "Read lengths are not consistent across profiles.\
                Make sure that variant profiles are generated with parameters\
                matching the profiles of samples on which they will be used."
                );
        if (!variantProfiles[i].variantStructureIsPresent())
        {
            if (params.variantFile == "")
            {
                std::string msg = "Variant description could not be loaded from profile. A variant description file needs to be specified.";
                throw std::runtime_error(msg.c_str());
            } else {
                if (!variantProfiles[i].loadVariantStructure(params.variantFile, variantProfiles[i].getName()))
                {
                    std::string msg =  "Description of variant " + variantProfiles[i].getName() + " could not be loaded from file " + params.variantFile + ".";
                    throw std::runtime_error(msg.c_str());
                } else {
                    variantProfiles[i].getFilter() = ReadPairFilter(variantProfiles[i].getVariant().getAllBreakpoints(), variantProfiles[i].getMargin(), 0);
                }
            }
        }
    }
    return;
}

inline void checkSampleParameters(std::vector<std::string> & sampleProfiles, int sMinVar, int sMaxVar, int readLengthVar, genotypeParameters & params)
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
        if (s.getLibraryDistribution().getReadLength() != readLengthVar)
            throw std::runtime_error("Read length in sample does not match read length in variant profiles.");
        if (s.getLibraryDistribution().getMinInsert() < sMinVar)
            throw std::runtime_error("Minimum insert size in sample distribution is not contained in at least one variant profile.");
        if (s.getLibraryDistribution().getMaxInsert() > sMaxVar)
            throw std::runtime_error("Maximum insert size in sample distribution is not contained in at least one variant profile.");
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
    if (! boost::filesystem::create_directory(directory)) {
        if (! boost::filesystem::is_directory(directory)) {
            std::cerr << "Could not create variant distribution directory" << std::endl;
            return;
        }
    }

    directory += (variantName + "/");
    if (! boost::filesystem::create_directory(directory)) {
        if (! boost::filesystem::is_directory(directory)) {
            std::cerr << "Could not create variant distribution directory" << std::endl;
            return;
        }
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

    return params;
}