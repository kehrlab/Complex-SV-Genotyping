#include "sample.hpp"
#include "bamFileHandler.hpp"
#include "custom_types.hpp"
#include "filter.hpp"
#include "genotypeResult.hpp"
#include "libraryDistribution.hpp"
#include "recordManager.hpp"
#include "variant.hpp"
#include "variantGenotyper.hpp"
#include <chrono>
#include <boost/filesystem.hpp>
#include <cstring>
#include <ostream>
#include <stdexcept>
#include <unordered_map>

// Sample::Sample(Sample && s): 
//     filename(s.filename), sampleName(s.sampleName),
//     minMapQ(s.minMapQ),
//     distributionDirectory(s.distributionDirectory), sampledRegions(s.sampledRegions), 
//     sampleDistribution(s.sampleDistribution)
// {
//     this->bamFileOpen = false;
// }

Sample & Sample::operator=(Sample s)
{
    this->filename = s.filename;
    this->sampleName = s.sampleName;
    this->minMapQ = s.minMapQ;
    this->distributionDirectory = s.distributionDirectory;
    this->sampledRegions = s.sampledRegions;
    this->sampleDistribution = s.sampleDistribution;
    this->bamFileOpen = false;
    
    return *this;
}

Sample::Sample(
    std::string filename,
    const std::vector<GenomicRegion> & givenRegions
) : filename(filename), sampledRegions(givenRegions)
{
    this->distributionDirectory = "";
    this->sampleName = this->filename;
    this->minMapQ = 0;
    this->bamFileOpen = false;

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
    this->chromosomeLengths = bamFile.getContigLengths();
    this->sampleName = bamFile.getSampleName();
}

// void Sample::createDistributionDirectory()
// {
//     if (!this->options.isOptionOutputDistributions())
//         return;

//     this->distributionDirectory = this->filename + "_distributions";

//     if (! boost::filesystem::create_directory(this->distributionDirectory))
//     {
//         if (! boost::filesystem::is_directory(this->distributionDirectory))
//         {
//             std::string msg = "Could not create directory '" + this->distributionDirectory + "', but flag -d was set. Make sure you have write permissions in the directory containing bam files.";
//             throw std::runtime_error(msg.c_str());
//         }
//     }
// }

void Sample::calculateDefaultDistributions()
{
    std::unordered_map<std::string, TemplatePosition> insertPositions;

    int readLength;
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

ContigInfo Sample::getContigInfo()
{
    return this->bamFile.getContigInfo();
}

std::unordered_map<std::string, int> Sample::getContigLengths()
{
    return this->bamFile.getContigLengths();
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

    // write version
    stream.write("0.9", 3);

    // write sample name
    int sampleNameLen = this->sampleName.size() + 1;
    stream.write(reinterpret_cast<const char *>(&sampleNameLen), sizeof sampleNameLen);
    stream.write(this->sampleName.c_str(), sampleName.size());
    stream.write("\0", 1);

    // write sample path
    int samplePathLength = this->filename.size() + 1;
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
    int nChrom = this->chromosomeLengths.size();
    stream.write(reinterpret_cast<const char *>(&nChrom), sizeof nChrom);
    for (auto & chr : this->chromosomeLengths)
    {
        int cNameLen = chr.first.size() + 1;
        stream.write(reinterpret_cast<const char *>(&cNameLen), sizeof cNameLen);
        stream.write(chr.first.c_str(), cNameLen - 1);
        stream.write("\0", 1);
        stream.write(reinterpret_cast<const char *>(&chr.second), sizeof chr.second);
    }

    // write regions used for distribution creation
    int nRegions = this->regionStrings.size();
    stream.write(reinterpret_cast<const char *>(&nRegions), sizeof nRegions);
    for (auto & r : this->regionStrings)
    {
        int rNameLen = r.size() + 1;
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

    char * buffer = new char[100];
    // magic string
    char * magic = new char[10];
    stream.read(magic, 10);
    
    if (!stream.good())
        throw std::runtime_error("Unable to read magic string.");
    if (std::strcmp(magic, "GENOTYPER\1") != 0)
        throw std::runtime_error("Magic string is wrong");
    
    // version
    char * version = new char[3];
    stream.read(version, 3);

    // sample name
    int sampleNameSize;
    stream.read(reinterpret_cast<char *>(&sampleNameSize), sizeof(int));
    char * tempString = new char[sampleNameSize];
    stream.read(tempString, sampleNameSize);
    if (tempString[sampleNameSize - 1] != '\0')
        throw std::runtime_error("In Sample Profile: String must end with 0!");
    this->sampleName = std::string(tempString);
    delete[] tempString;
    
    // sample path
    int samplePathSize;
    stream.read(reinterpret_cast<char *>(&samplePathSize), sizeof(int));
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
    std::unordered_map<std::string, int>().swap(this->chromosomeLengths);
    int nChrom;
    stream.read(reinterpret_cast<char * >(&nChrom), sizeof(int));
    for (int i = 0; i < nChrom; ++i)
    {
        int l;
        int cLength;

        stream.read(reinterpret_cast<char *>(&l), sizeof(int));

        char * tempString = new char[l];
        stream.read(reinterpret_cast<char *>(tempString), l);
        std::string cName = std::string(tempString);
        delete[] tempString;
        
        stream.read(reinterpret_cast<char *>(&cLength), sizeof(int));
        this->chromosomeLengths[cName] = cLength;
    }

    // regions
    std::vector<std::string>().swap(this->regionStrings);
    int nRegions;
    stream.read(reinterpret_cast<char * >(&nRegions), sizeof(int));
    for (int i = 0; i < nRegions; ++i)
    {
        int l;
        stream.read(reinterpret_cast<char *>(&l), sizeof(int));

        char * tempString = new char[l];
        stream.read(tempString, l);
        std::string region = std::string(tempString);
        delete[] tempString;
        this->regionStrings.push_back(region);
    }

    // insert size values
    std::vector<float> distribution(this->sampleDistribution.getMaxInsert() - this->sampleDistribution.getMinInsert() + 1);
    for (int i = 0; i < distribution.size(); ++i)
        stream.read(reinterpret_cast<char *>(&distribution[i]), sizeof(float));
    this->sampleDistribution.getInsertDistribution() = distribution;

    delete[] magic;
    delete[] version;
    delete[] buffer;

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
    for (auto & chr : this->chromosomeLengths)
        std::cout << "\t" << chr.first << "\t" << chr.second << std::endl;
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