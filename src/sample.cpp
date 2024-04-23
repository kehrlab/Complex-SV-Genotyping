#include "sample.hpp"

Sample::Sample(
    std::string filename,
    const std::vector<GenomicRegion> & givenRegions
) : filename(filename), sampledRegions(givenRegions)
{
    this->distributionDirectory = "";
    this->sampleName = this->filename;
    this->minMapQ = 0;
    this->bamFileOpen = false;

    this->regionStrings = std::vector<std::string>();
    for (auto & r : this->sampledRegions)
	    this->regionStrings.push_back(r.getRegionString());

    openBamFile();
    calculateDefaultDistributions();
    closeBamFile();
}

Sample::Sample(std::string profilePath)
{
    this->distributionDirectory = "";
    this->bamFileOpen = false;

    std::ifstream stream;
    this->sampleDistribution = LibraryDistribution();
    readSampleProfile(profilePath);
}

void Sample::openBamFile()
{
    this->bamFile.open(this->filename);
    this->bamFileOpen = true;
    this->contigInfo = bamFile.getContigInfo();

    this->sampleName = this->bamFile.getSampleName();
    if (this->sampleName == "")
    {
        std::string file = this->filename.substr(this->filename.find_last_of("/") + 1);
        this->sampleName = file.substr(0, file.find(".bam"));
    }
}

void Sample::calculateDefaultDistributions()
{
    std::unordered_map<std::string, TemplatePosition> insertPositions;

    int readLength = 0;
    if (this->sampledRegions.size() == 0) 
        insertPositions = this->bamFile.get_insert_size_positions(readLength);
    else 
        insertPositions = this->bamFile.get_insert_size_positions(readLength, this->regionStrings);

    this->sampleDistribution = LibraryDistribution(insertPositions);
    this->sampleDistribution.getReadLength() = readLength;

    std::unordered_map<std::string, TemplatePosition>().swap(insertPositions);
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

const ContigInfo & Sample::getContigInfo()
{
    return this->contigInfo;
}

std::unordered_map<std::string, int32_t> Sample::getContigLengths()
{
    return this->contigInfo.getContigLengths();
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

int Sample::getMaxReadLength()
{
    return this->sampleDistribution.getReadLength();
}

void Sample::writeSampleProfile(std::string profilePath)
{
    std::ofstream stream(profilePath, std::ios_base::binary | std::ios_base::out);
    if (!stream.is_open())
    {
        std::string msg = "Could not open profile path (" + profilePath + ") for reading.";
        throw std::runtime_error(msg.c_str());
    }
    
    // magic string
    stream.write("GENOTYPER\1", 10);

    // write sample name
    uint32_t sampleNameLen = this->sampleName.size() + 1;
    stream.write(reinterpret_cast<const char *>(&sampleNameLen), sizeof sampleNameLen);
    stream.write(this->sampleName.c_str(), sampleName.size());
    stream.write("\0", 1);

    // write sample path
    uint32_t samplePathLength = this->filename.size() + 1;
    stream.write(reinterpret_cast<const char *>(&samplePathLength), sizeof samplePathLength);
    stream.write(this->filename.c_str(), this->filename.size());
    stream.write("\0", 1);

    // write read length and number of read pairs considered
    stream.write(reinterpret_cast<const char *>(& this->sampleDistribution.getReadLength()), sizeof(int));
    stream.write(reinterpret_cast<const char *>(& this->sampleDistribution.getNumReadPairs()), sizeof(int));

    // write min mapQ
    stream.write(reinterpret_cast<const char *>(& this->minMapQ), sizeof this->minMapQ);

    // write sMin and sMax
    stream.write(reinterpret_cast<const char *>(&this->sampleDistribution.getMinInsert()), sizeof(int));
    stream.write(reinterpret_cast<const char *>(&this->sampleDistribution.getMaxInsert()), sizeof(int));

    // write MEAN and SD
    stream.write(reinterpret_cast<const char *>(&this->sampleDistribution.getInsertMean()), sizeof(float));
    stream.write(reinterpret_cast<const char *>(&this->sampleDistribution.getInsertSD()), sizeof(float));

    // write chromosome names and lengths
    // write number of chromosomes
    uint32_t nChrom = this->contigInfo.cNames.size();
    stream.write(reinterpret_cast<const char *>(&nChrom), sizeof nChrom);
    for (uint32_t i = 0; i < nChrom; ++i)
    {
        uint32_t cNameLen = this->contigInfo.cNames[i].size() + 1;
        stream.write(reinterpret_cast<const char *>(&cNameLen), sizeof cNameLen);
        stream.write(this->contigInfo.cNames[i].c_str(), cNameLen - 1);
        stream.write("\0", 1);
        stream.write(reinterpret_cast<const char *>(&this->contigInfo.cLengths[i]), sizeof(int32_t));
    }

    // write regions used for distribution creation
    uint32_t nRegions = this->regionStrings.size();
    stream.write(reinterpret_cast<const char *>(&nRegions), sizeof nRegions);
    for (auto & r : this->regionStrings)
    {
        uint32_t rNameLen = r.size() + 1;
        stream.write(reinterpret_cast<const char *>(&rNameLen), sizeof rNameLen);
        stream.write(r.c_str(), rNameLen - 1);
        stream.write("\0", 1);
    }

    // write actual distribution
    for (auto & p : this->sampleDistribution.getInsertDistribution())
        stream.write(reinterpret_cast<const char *>(&p), sizeof(float));
    
    stream.close();
}

void Sample::readSampleProfile(std::string profilePath)
{
    std::ifstream stream(profilePath, std::ios_base::binary | std::ios_base::out);
    if (!stream.is_open())
    {
        std::string msg = "Could not open profile path (" + profilePath + ") for reading.";
        throw std::runtime_error(msg.c_str());
    }
    
    // compare magic string
    char * tempString = new char[10];
    stream.read(tempString, 10);
    if (!stream.good())
        throw std::runtime_error("Unable to read magic string.");
    if (tempString[9] != '\1')
        throw std::runtime_error("Magic string does not match. Wrong / corrupted profile?");
    tempString[9] = '\0';
    if (std::strcmp(tempString, "GENOTYPER") != 0)
        throw std::runtime_error("Magic string does not match. Wrong / corrupted profile?");
    delete[] tempString;

    // sample name
    uint32_t sampleNameSize;
    stream.read(reinterpret_cast<char *>(&sampleNameSize), sizeof(uint32_t));
    tempString = new char[sampleNameSize];
    stream.read(tempString, sampleNameSize);
    if (tempString[sampleNameSize - 1] != '\0')
        throw std::runtime_error("In Sample Profile: String must end with 0!");
    this->sampleName = std::string(tempString);
    delete[] tempString;
    
    // sample path
    uint32_t samplePathSize;
    stream.read(reinterpret_cast<char *>(&samplePathSize), sizeof(uint32_t));
    tempString = new char[samplePathSize];
    stream.read(tempString, samplePathSize);
    if (tempString[samplePathSize - 1] != '\0')
        throw std::runtime_error("In Sample Profile: String must end with 0!");
    this->filename = std::string(tempString);
    delete[] tempString;

    // max read length and number of reads
    stream.read(reinterpret_cast<char *>(&this->sampleDistribution.getReadLength()), sizeof(int));
    stream.read(reinterpret_cast<char *>(&this->sampleDistribution.getNumReadPairs()), sizeof(int));
    // minMapQ and insert size numbers
    stream.read(reinterpret_cast<char *>(&this->minMapQ), sizeof(int));
    stream.read(reinterpret_cast<char *>(&this->sampleDistribution.getMinInsert()), sizeof(int));
    stream.read(reinterpret_cast<char *>(&this->sampleDistribution.getMaxInsert()), sizeof(int));
    stream.read(reinterpret_cast<char *>(&this->sampleDistribution.getInsertMean()), sizeof(float));
    stream.read(reinterpret_cast<char *>(&this->sampleDistribution.getInsertSD()), sizeof(float));

    // chromosome names and lengths
    this->contigInfo = ContigInfo();
    uint32_t nChrom;
    stream.read(reinterpret_cast<char * >(&nChrom), sizeof(uint32_t));
    for (uint32_t i = 0; i < nChrom; ++i)
    {
        uint32_t l;
        int cLength;

        stream.read(reinterpret_cast<char *>(&l), sizeof(uint32_t));

        char * tempString = new char[l];
        stream.read(reinterpret_cast<char *>(tempString), l);
        std::string cName = std::string(tempString);
        delete[] tempString;
        
        stream.read(reinterpret_cast<char *>(&cLength), sizeof(int32_t));
        this->contigInfo.cNames.push_back(cName);
        this->contigInfo.cLengths.push_back(cLength);
    }
    this->contigInfo.calculateGlobalContigPositions();

    // regions
    std::vector<std::string>().swap(this->regionStrings);
    uint32_t nRegions;
    stream.read(reinterpret_cast<char * >(&nRegions), sizeof(uint32_t));
    for (uint32_t i = 0; i < nRegions; ++i)
    {
        uint32_t l;
        stream.read(reinterpret_cast<char *>(&l), sizeof(uint32_t));

        char * tempString = new char[l];
        stream.read(tempString, l);
        std::string region = std::string(tempString);
        delete[] tempString;
        this->regionStrings.push_back(region);
    }

    // insert size values
    std::vector<float> distribution(this->sampleDistribution.getMaxInsert() - this->sampleDistribution.getMinInsert() + 1);
    for (uint32_t i = 0; i < distribution.size(); ++i)
        stream.read(reinterpret_cast<char *>(&distribution[i]), sizeof(float));
    this->sampleDistribution.getInsertDistribution() = distribution;

    stream.close();
}

void Sample::printSampleProfile()
{
    std::cout << "--------------------------------------------------------------" << std::endl;
    std::cout << "Version: 0.9" << std::endl;
    std::cout << "Sample Name: " << this->sampleName << std::endl;
    std::cout << "Filepath: " << this->filename << std::endl;
    std::cout << "Read Length: " << this->sampleDistribution.getReadLength() << std::endl;
    std::cout << "Read Number: " << this->sampleDistribution.getNumReadPairs() << std::endl;
    std::cout << "Min MapQ: " << this->minMapQ << std::endl;
    std::cout << "sMin: " << this->sampleDistribution.getMinInsert() << std::endl;
    std::cout << "sMax: " << this->sampleDistribution.getMaxInsert() << std::endl;
    std::cout << "Mean: " << this->sampleDistribution.getInsertMean() << std::endl;
    std::cout << "SD: " << this->sampleDistribution.getInsertSD() << std::endl;
    std::cout << "Chromosomes: "<< std::endl;
    for (uint32_t i = 0; i < this->contigInfo.cNames.size(); ++i)
        std::cout << "\t" << this->contigInfo.cNames[i] << "\t" << this->contigInfo.cLengths[i] << std::endl;
    std::cout << "Regions:" << std::endl;
    for (auto & r : this->regionStrings)
        std::cout << "\t" << r << std::endl;
    std::cout << "Distribution: " << std::endl;
    for (auto & p : this->sampleDistribution.getInsertDistribution())
        std::cout << p << "\t";
    std::cout << std::endl;
    std::cout << "--------------------------------------------------------------" << std::endl;
}

std::string Sample::getSampleName()
{
    return this->sampleName;
}
