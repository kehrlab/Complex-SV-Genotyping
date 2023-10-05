#ifndef ALLELEHEADER
#define ALLELEHEADER

#include "breakpoint.hpp"
#include "junction.hpp"
#include "seqFileHandler.hpp"
#include "vMapManager.hpp"
#include <unordered_set>

class Allele
{
    std::vector<Junction> novelJunctions;
    std::vector<Breakpoint> referenceBreakpoints;
    std::string alleleName;

    std::vector<std::string> chromosomeNames;
    std::vector<VariantMapManager> chromosomeMaps;

    seqan::String<seqan::Dna5String> chromosomeSequences;
    seqan::String<seqan::CharString> sequenceIDs;

    void convertJunctionsToBreakpoints();

    void splitAllJunctionsToBreakpoints();
    void removeDuplicateBreakpoints();
    void reorderBreakPoints();
    void joinBreakpoints();
    std::unordered_set<std::string> getUniqueBpRefNames();
    void sortBreakpointsOnReference(std::string);
    void insertBreakpoint(Breakpoint & bp, JunctionRegion & jRegion);

    public:
    Allele();
    Allele(std::vector<Breakpoint>);
    Allele(std::vector<Junction> novelJunctions);
    Allele(std::vector<Junction> novelJunctions, std::string);

    std::vector<Junction> & getNovelJunctions();
    std::vector<Breakpoint> & getReferenceBreakpoints();
    std::string getName();
    seqan::String<seqan::CharString> getSequenceIDs();
    seqan::String<seqan::Dna5String> getChromosomeSequences();
    void print();
    void createChromosomeMaps(std::vector<Breakpoint> &, int, std::unordered_map<std::string, int>);
    void createAlleleSequences(SeqFileHandler &);

    std::unordered_map<std::string, JunctionRegion> getChromosomeStructures(std::vector<Breakpoint> &);
    
    std::vector<Junction> getJunctionsOnChromosome(std::string);
    void createChromosomeMap(std::string, std::vector<Breakpoint> &, int, std::unordered_map<std::string, int>);

    void writeSequence(SeqFileHandler &, std::string);

    std::vector<VariantMapManager> & getChromosomeMaps();
    std::vector<std::string> & getChromosomeNames();
    VariantMapManager & getChromosomeMap(std::string);
    void clearSequences();
};

#endif
