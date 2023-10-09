#include "genotyper.hpp"
#include "genotypeResult.hpp"
#include "libraryDistribution.hpp"
#include "regionSampler.hpp"
#include "variantProfile.hpp"
#include <chrono>
#include <ios>
#include <unordered_map>
#include <unordered_set>

Genotyper::Genotyper()
{
    this->maxFilterMargin = 500;
}

Genotyper::Genotyper(ProgramOptions options)
{
    this->maxFilterMargin = 500;
    setOptions(options);
    getBamFileNames();
    createVariants();
    extractContigInfos();
    sampleRegionsForInsertSizeDistribution();
    this->options.determineNumThreads(this->fileNames.size(), this->variants.size());
    createSamples();

    if (this->options.isOptionCreateVariantProfiles())
        createVariantProfiles();
    else
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
}

void Genotyper::createVariants()
{
    variantParser vParser(this->options.getVariantFileName());
    std::vector<variantData> allVariantJunctions = vParser.getVariantJunctions();
    std::vector<std::string> variantNames = vParser.getVariantNames();
    std::vector<std::vector<std::string>> alleleNames = vParser.getAlleleNames();

    for (int i = 0; i < variantNames.size(); ++i)
        this->variants.push_back(complexVariant(variantNames[i], alleleNames[i], allVariantJunctions[i], this->options.getVariantFileName()));
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
    openReferenceFile();
    std::vector<GenomicRegion> regionsForSampling = this->options.getSamplingRegions();
    if (this->contigInfos.size() > 0 && !this->options.isOptionUseWholeGenome())
        this->regionSampler = RegionSampler(this->contigInfos, this->referenceFile, this->fileNames, this->options, regionsForSampling);
    closeReferenceFile();
}

void Genotyper::createSamples()
{   
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    
    sampleDistributions = std::vector<LibraryDistribution>(this->fileNames.size());
    std::vector<int> readLengths(this->fileNames.size(), 0);

    // open samples and calculate library distributions
    for (int i = 0; i < this->fileNames.size(); ++i) 
    {
        Sample s(this->fileNames[i], this->options.getRefFileName(), this->options, this->regionSampler);
        readLengths[i] = s.getMaxReadLength();
        this->sampleDistributions[i] = s.getLibraryDistribution();
    }

    // determine maximum read length across all samples
    readLengths.size() > 0 ? this->maxReadLength = 0 : this->maxReadLength = 150;
    for (int l : readLengths)
        this->maxReadLength = std::max(this->maxReadLength, l);

    // reset the region sampler, as it is no longer needed
    this->regionSampler = RegionSampler();

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    if (this->options.isOptionProfile())
        std::cout << "Total time taken to create default distributions for all files: " << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << "s" << std::endl;
}

void Genotyper::createVariantProfiles()
{
    int overlap = 20;
    int sMin {1}, sMax {1000};

    for (int i = 0; i < this->variants.size(); ++i) 
        this->variantProfiles.push_back(VariantProfile(this->variants[i], this->maxFilterMargin, overlap, this->maxReadLength, sMin, sMax, this->contigLengths, this->options));
    
    #pragma omp parallel for num_threads(this->options.getNumberOfThreads())
    for (int i = 0; i < this->variantProfiles.size(); ++i) 
       this->variantProfiles[i].calculateAlleleMasks();

    // this->variantProfiles[0].writeProfile(std::string("variantTest.profile"));
    // this->variantProfiles[0].readProfile("variantTest.profile");
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
    else
    {
        openReferenceFile();
        if (this->referenceFile.isOpen())
            this->contigLengths = this->referenceFile.getSequenceLengths();
        closeReferenceFile();
    }
}

void Genotyper::createAlleleMaps()
{
    for (complexVariant & cSV : this->variants) 
        cSV.createAlleleMaps(this->maxFilterMargin, this->contigLengths);
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
    openReferenceFile();
    if (!this->referenceFile.isOpen())
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
    closeReferenceFile();
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
            if (this->variantProfiles.size() > 0)
                vG.genotype(this->variantProfiles[i]);
            else
                vG.genotype();
            this->genotypeResults[i][j] = vG.getResult();
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            if (this->options.isOptionProfile())
                std::cout << "Total time to genotype variant " << this->variants[i].getName() << " in file " << this->fileNames[j] << ": " << std::chrono::duration_cast<std::chrono::seconds>(end-begin).count() << "s" << std::endl;
        }
    }
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
