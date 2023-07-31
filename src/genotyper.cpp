#include "genotyper.hpp"
#include "genotypeResult.hpp"
#include "libraryDistribution.hpp"
#include "regionSampler.hpp"
#include <chrono>
#include <unordered_map>
#include <unordered_set>

Genotyper::Genotyper()
{}

Genotyper::Genotyper(ProgramOptions options)
{
    setOptions(options);
    getBamFileNames();
    openReferenceFile();
    createVariants();
    extractContigInfos();
    sampleRegionsForInsertSizeDistribution();
    this->options.determineNumThreads(this->fileNames.size(), this->variants.size());
    createSamples();
    createAlleleMaps();
    printVariantDetails();
    writeVariantAlleles();
}

void Genotyper::setOptions(ProgramOptions & options)
{
    this->options = options;
}

void Genotyper::getBamFileNames()
{
    this->fileNames = this->options.getBamFileNames();
}

void Genotyper::openReferenceFile()
{
    if (!this->referenceFile.isOpen() && this->options.getRefFileName() != "")
        this->referenceFile.open(this->options.getRefFileName());
    if (this->referenceFile.isOpen() && this->options.isOptionLoadToMemory() && !this->options.isOptionUseWholeGenome()) {
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        this->referenceFile.loadChromosomeSequences();
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        if (this->options.isOptionProfile())
            std::cout << "Time taken to load reference to memory: " << std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << "ms" << std::endl;
    }
}

void Genotyper::createVariants()
{
    variantParser vParser(this->options.getVariantFileName());
    std::vector<variantData> allVariantJunctions = vParser.getVariantJunctions();
    std::vector<std::string> variantNames = vParser.getVariantNames();
    std::vector<std::vector<std::string>> alleleNames = vParser.getAlleleNames();

    for (int i = 0; i < variantNames.size(); ++i)
        this->variants.push_back(complexVariant(variantNames[i], alleleNames[i], allVariantJunctions[i]));
}

void Genotyper::extractContigInfos()
{
    for (auto filename : this->fileNames)
    {
        BamFileHandler tempFileHandler(filename);
        this->contigInfos.push_back(tempFileHandler.getContigInfo());
        tempFileHandler.closeInputFile();
    }
    gatherContigLengths();
}

void Genotyper::sampleRegionsForInsertSizeDistribution()
{
    std::vector<GenomicRegion> regionsForSampling = this->options.getSamplingRegions();
    if (this->contigInfos.size() > 0 && !this->options.isOptionUseWholeGenome())
        this->regionSampler = RegionSampler(this->contigInfos, this->referenceFile, this->fileNames, this->options, regionsForSampling);
}

void Genotyper::createSamples()
{   
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    for (int i = 0; i < this->fileNames.size(); ++i)
        sampleDistributions.push_back(LibraryDistribution());
    std::vector<int> margins(this->fileNames.size(), 0);

    //#pragma omp parallel for num_threads(this->options.getNumberOfThreads())
    for (int i = 0; i < this->fileNames.size(); ++i) 
    {
        Sample s(this->fileNames[i], this->referenceFile, this->options, this->regionSampler, this->variants);
        s.open();
        margins[i] = s.getFilterMargin();
        this->sampleDistributions[i] = s.getLibraryDistribution();
	s.close();
    }
    margins.size() > 0 ? this->maxFilterMargin = 0 : this->maxFilterMargin = 1000;
    for (int m : margins)
    	this->maxFilterMargin = std::max(this->maxFilterMargin, m);

    // calculateMaxFilterMargin(samples); // could make this one sample-specific
    this->regionSampler = RegionSampler();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    if (this->options.isOptionProfile())
        std::cout << "Total time taken to create default distributions for all files: " << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << "s" << std::endl;
}

void Genotyper::calculateMaxFilterMargin(std::vector<Sample> & samples)
{
    samples.size() > 0 ? this->maxFilterMargin = 0 : this->maxFilterMargin = 1000;
    for (Sample & s : samples)
        this->maxFilterMargin = std::max(this->maxFilterMargin, s.getFilterMargin());
}

void Genotyper::gatherContigLengths()
{
    if (this->contigInfos.size() > 0)
    {
        for (ContigInfo & cInfo : contigInfos)
        {
            for (int i = 0; i < cInfo.cNames.size(); ++i)
            {
                if (this->contigLengths.find(cInfo.cNames[i]) != this->contigLengths.end())
                    this->contigLengths[cInfo.cNames[i]] = std::min(this->contigLengths[cInfo.cNames[i]], cInfo.cLengths[i]);
                else
                    this->contigLengths[cInfo.cNames[i]] = cInfo.cLengths[i];
            }
        }
    }
    else if (this->referenceFile.isOpen())
    {
        this->contigLengths = this->referenceFile.getSequenceLengths();
    }
}

void Genotyper::createAlleleMaps()
{
    for (complexVariant & cSV : this->variants) 
        for (Allele & allele : cSV.getAlleles()) 
            allele.createChromosomeMaps(cSV.getAllBreakpoints(), this->maxFilterMargin, this->contigLengths, this->referenceFile); 
}

void Genotyper::printVariantDetails()
{
    if (this->options.isOptionVerbose())
    {
        std::cout << std::endl << "Specified Variants: " << std::endl << std::endl;
        for (auto& cSV : this->variants) {
            cSV.print();
            std::cout << std::endl << std::endl;
        }
    }
}

void Genotyper::writeVariantAlleles()
{
    if (!this->referenceFile.isOpen() || this->options.getRefFileName() == "")
        return;
    if (this->options.getSequenceDirectory() == "")
        return;

    if (! boost::filesystem::create_directory(this->options.getSequenceDirectory()))
    {
        if (! boost::filesystem::is_directory(this->options.getSequenceDirectory()))
        {
            std::cout << "Could not create directory '" + this->options.getSequenceDirectory() + "'." << std::endl;
            return;
        }
    }
    for (complexVariant & variant : this->variants)
    {
        std::string variantDirName = this->options.getSequenceDirectory() + "/" + variant.getName();
        if (! boost::filesystem::create_directory(variantDirName))
        {
            if (!boost::filesystem::is_directory(variantDirName)) {
                std::cout << "Could not create directory '" + variantDirName + "'." << std::endl;
                return;
            }
        }
        for (Allele & allele : variant.getAlleles()) {
            allele.createAlleleSequences(this->referenceFile);
            allele.writeSequence(this->referenceFile, variantDirName);
            allele.clearSequences();
        }
    }
}

void Genotyper::genotypeAllSamples()
{	
    createResultVector();
    
    #pragma omp parallel for schedule(dynamic) collapse(2) num_threads(this->options.getNumberOfThreads())
    for (int j = 0; j < this->fileNames.size(); ++j)
    {
        for (int i = 0; i < this->variants.size(); ++i)
        {
            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            VariantGenotyper vG(this->variants[i], this->fileNames[j], this->sampleDistributions[j], this->options);
            vG.genotype();
            this->genotypeResults[i][j] = vG.getResult();
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            if (this->options.isOptionProfile())
                std::cout << "Total time to genotype variant " << this->variants[i].getName() << " in file " << this->fileNames[j] << ": " << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << "s" << std::endl;
        }
    }
    //closeSamples();
    closeReferenceFile();
}

void Genotyper::createResultVector()
{
    for (int i = 0; i < this->variants.size(); ++i)
    {
        std::vector<GenotypeResult> variantResults(this->fileNames.size());
        this->genotypeResults.push_back(variantResults);
    }
}

void Genotyper::writeToVCF()
{
    if (this->options.getVCFName() != "")
    {
        VcfWriter vcfWriter(this->options.getVCFName(), this->variants, this->fileNames, this->genotypeResults);
        vcfWriter.write();
    }
}

void Genotyper::writeResults()
{
    std::string allFileName = this->options.getOutputFileName() + "_allVariants.tsv";
    std::ofstream fAll;

    if (this->options.isOptionWriteToFile()) 
    {
        fAll.open(allFileName);
        if (fAll.is_open())
            fAll << "Variant\tFile\tGenotype\tMean Quality\tLower Bound\tUpper Bound\tCertainty\t#Reads\tAvgMapQ\tMinMapQ\tMaxMapQ\tQualityPass" << std::endl;
    }
    
    for (int i = 0; i < this->variants.size(); ++i)
    {
        std::string filename = this->options.getOutputFileName() + "_" + this->variants[i].getName() + ".tsv";
        if (this->genotypeResults[i].size() == 0)
            continue;
        
        if (this->options.isOptionWriteToFile())
        {
            std::ofstream f(filename);

            if (f.is_open())
            {
                f << "File\tGenotype\tMean Quality\tLower Bound\tUpper Bound\tCertainty\t#Reads\tAvgMapQ\tMinMapQ\tMaxMapQ\tQualityPass" << std::endl;
                for (auto result : this->genotypeResults[i])
                    f << result.getOutputString();
            }
            f.close();
            if (fAll.is_open())
            {
                for (auto result : this->genotypeResults[i]) {
                    fAll << this->variants[i].getName() << "\t";
                    fAll << result.getOutputString();
                }
            }
        } else {
            std::cout << this->variants[i].getName() << std::endl;
            std::cout << "File\tGenotype\tMean Quality\tLower Bound\tUpper Bound\tCertainty\t#Reads\tAvgMapQ\tMinMapQ\tMaxMapQ\tQualityPass" << std::endl;
            for (auto result : this->genotypeResults[i])
                std::cout << result.getOutputString();
            std::cout << std::endl;
        }
    }
    if (fAll.is_open())
        fAll.close();
}

void Genotyper::writeStats()
{
    if (!this->options.isOptionStats())
        return;

    std::string statsFileName = this->options.getOutputFileName() + "_stats.tsv";
    std::ofstream fStats(statsFileName);

    if (!fStats.is_open())
        return;

    std::unordered_set<std::string> groups;
    std::vector<std::unordered_map<std::string, float>> readPairGroupQualities;

    for (int i = 0; i < this->variants.size(); ++i)
    {
        if (this->genotypeResults[i].size() == 0)
            continue;
        for (GenotypeResult & result : this->genotypeResults[i])
        {
            std::unordered_map<std::string, float> groupQualities = result.getReadStats();
            for (auto key : groupQualities)
                groups.insert(key.first);
            readPairGroupQualities.push_back(groupQualities);
        }
    }

    fStats << "File\tVariant\t";
    for (auto s : groups)
        fStats << s << "\t";
    fStats << std::endl;

    int index = 0;
    for (int i = 0; i < this->variants.size(); ++i)
    {
        if (this->genotypeResults[i].size() == 0)
            continue;
        for (GenotypeResult & result : this->genotypeResults[i])
        {
            fStats << result.getFilename() << "\t" << this->variants[i].getName() << "\t";
            for (auto s : groups)
            {
                if (readPairGroupQualities[index].find(s) != readPairGroupQualities[index].end())
                    fStats << readPairGroupQualities[index][s];
                fStats << "\t";
            }
            fStats << std::endl;
            ++index;
        }
    }
    fStats.close();
}

void Genotyper::closeReferenceFile()
{
    this->referenceFile.close();
}
