#include "seqFileHandler.hpp"

SeqFileHandler::SeqFileHandler()
{
    this->filename = "";
    this->success = false;
    this->sequencesInMemory = false;
}

SeqFileHandler::SeqFileHandler(std::string filename)
{
    this->sequencesInMemory = false;
    setFileName(filename);
    createIndex();
    if (!this->success)
        std::cout << "Could not build index!" << std::endl;
    if (!openFile())
        throw std::runtime_error("Could not open sequence file / index.");
}

void SeqFileHandler::setFileName(std::string filename)
{
    this->filename = filename;
}

void SeqFileHandler::createIndex()
{
    std::string faiFileName = this->filename + ".fai";
    this->success = seqan::build(this->faiIndex, filename.c_str(), faiFileName.c_str());
    seqan::save(faiIndex);
}

bool SeqFileHandler::openFile()
{
    return seqan::open(this->faiIndex, this->filename.c_str());
}

int SeqFileHandler::getIndex(std::string rName)
{
    int idx = -1;
    if (!seqan::getIdByName(idx, this->faiIndex, rName)) {
        idx = -1;
        std::cout << "ERROR: FaiIndex has no entry for " + rName << std::endl;
    }
    return idx;
}

seqan::Dna5String SeqFileHandler::getBaseAtPosition(std::string rName, int position)
{
    seqan::Dna5String sequence;
    if (position >= getChrLength(rName))
        return sequence;

    if (this->sequencesInMemory && this->chromosomeSequences.find(rName) != this->chromosomeSequences.end()) {
        if (position < seqan::length(this->chromosomeSequences[rName][position]))
            sequence = this->chromosomeSequences[rName][position];
    }
    else
    {
        int idx = getIndex(rName);
        if (idx < 0)
            return sequence;
        seqan::readRegion(sequence, this->faiIndex, idx, position + 1, position + 2);
    }
    return sequence;
}

seqan::Dna5String SeqFileHandler::getRegionSequence(std::string rName, int begin, int end)
{
    seqan::Dna5String sequence;
    if (this->sequencesInMemory && this->chromosomeSequences.find(rName) != this->chromosomeSequences.end())
    {
        for (int i = begin; i <= end; ++i)
            if (i < seqan::length(this->chromosomeSequences[rName]))
                seqan::append(sequence, this->chromosomeSequences[rName][i]);
            else
                break;
    } 
    else
    {
        int idx = getIndex(rName);
        if (idx < 0)
            return sequence;
        ++begin;
        ++end;
        end = std::min(end + 1, getChrLength(rName));
        seqan::readRegion(sequence, this->faiIndex, idx, begin, end);
    }
    return sequence;
}

seqan::Dna5String SeqFileHandler::getChromosomeSequence(std::string rName)
{
    seqan::Dna5String sequence;
    
    if (this->sequencesInMemory && this->chromosomeSequences.find(rName) != this->chromosomeSequences.end())
    {
        return this->chromosomeSequences[rName];
    }
    else
    {
        int idx = getIndex(rName);
        if (idx < 0)
            return sequence;
        seqan::readSequence(sequence, this->faiIndex, idx);
    }
    return sequence;
}

void SeqFileHandler::writeSequencesToFile(seqan::String<seqan::CharString>ids, seqan::String<seqan::Dna5String> sequences, std::string outFileName)
{
    seqan::SeqFileOut outputFile(outFileName.c_str());
    if (seqan::isOpen(outputFile)) {
        try {
            #pragma omp critical
            writeRecords(outputFile, ids, sequences);
        }
        catch (seqan::Exception const & e)
        {
            std::cout << "ERROR: " << e.what() << std::endl;
        }
        seqan::close(outputFile);
    } else {
        std::cerr << "Could not open output file " << outFileName << std::endl;
    }
}

bool SeqFileHandler::containsChrSequence(std::string cName)
{
    int idx = getIndex(cName);
    if (idx == -1)
        return false;
    return true;
}

int SeqFileHandler::getChrLength(std::string cName)
{
    int idx = getIndex(cName);
    if (idx != -1)
        return seqan::sequenceLength(this->faiIndex, idx);
    return 0;
}

void SeqFileHandler::open(std::string filename)
{
    setFileName(filename);
    createIndex();
    if (!this->success)
        std::cout << "Could not build index!" << std::endl;
    if (!openFile()) {
        this->success = false;
        throw std::runtime_error("Could not open sequence file / index.");
    } else {
        this->success = true;
    }
}

bool SeqFileHandler::isOpen()
{
    return this->success;
}

void SeqFileHandler::close()
{
    if (seqan::isOpen(this->seqFileIn))
        seqan::close(this->seqFileIn);
}

std::unordered_map<std::string, int> SeqFileHandler::getSequenceLengths()
{
    std::unordered_map<std::string, int> seqLens;
    if (isOpen())
        for (unsigned i = 0; i < seqan::numSeqs(this->faiIndex); ++i)
            seqLens[std::string(seqan::toCString(seqan::sequenceName(this->faiIndex, i)))] = seqan::sequenceLength(this->faiIndex, i);
    return seqLens;
}

void SeqFileHandler::loadChromosomeSequences()
{
    for (auto it : getSequenceLengths())
        if (it.first.size() < 6)
            this->chromosomeSequences[it.first] = getChromosomeSequence(it.first);
    this->sequencesInMemory = true;
}

bool SeqFileHandler::hasSequencesInMemory()
{
    return this->sequencesInMemory;
}