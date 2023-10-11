#ifndef SEQFILEHANDLERHEADER
#define SEQFILEHANDLERHEADER

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <stdexcept>
#include <unordered_set>

class SeqFileHandler
{
    std::string filename;
    seqan::SeqFileIn seqFileIn;
    seqan::FaiIndex faiIndex;
    std::unordered_map<std::string, seqan::Dna5String> chromosomeSequences;
    bool success;
    bool sequencesInMemory;

    void createIndex();
    void setFileName(std::string);
    bool openFile();

    public:
    SeqFileHandler();
    SeqFileHandler(std::string);

    int getIndex(std::string);
    seqan::Dna5String getBaseAtPosition(std::string, int);
    seqan::Dna5String getRegionSequence(std::string, int, int);
    seqan::Dna5String getChromosomeSequence(std::string);
    void writeSequencesToFile(seqan::String<seqan::CharString>, seqan::String<seqan::Dna5String>, std::string);
    bool containsChrSequence(std::string);
    int getChrLength(std::string);
    void open(std::string);
    bool isOpen();
    void close();
    void loadChromosomeSequences();
    bool hasSequencesInMemory();
    std::unordered_map<std::string, int> getSequenceLengths();
};

#endif