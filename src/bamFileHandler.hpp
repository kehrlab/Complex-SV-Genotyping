#ifndef BAMFILEHANDLERHEADER
#define BAMFILEHANDLERHEADER

#include <iostream>
#include <seqan/bam_io.h>
#include "genomicRegion.hpp"
#include <unordered_map>
#include <vector>
#include "custom_types.hpp"
#include "options.hpp"
#include <htslib/sam.h>
#include "record.hpp"

class BamFileHandler
{
    std::string bamFileName;
    std::string baiFileName;
    std::vector<std::string> regions;
    int beginPos, endPos, rID;
    bool hasAlignments;
    int minQ = 0;

    samFile * bamFileIn;
    hts_idx_t *index;
    bam_hdr_t *header;

    ContigInfo contigInfo;

    public:
    BamFileHandler();
    BamFileHandler(std::string);
    BamFileHandler(std::string, int);
    BamFileHandler(std::string, ProgramOptions &);
    
    void open(std::string, ProgramOptions &);
    void open(std::string);
    void openInputFiles();
    void inferBaiFileName();
    void openBamFile();
    void openBaiFile();
    void readHeader();
    void closeInputFile();
    void setRegions(std::vector<GenomicRegion>);
    void resetRegions();
    std::string getFileName();
    void extractContigInfo();

    // get_insert_size_positions
    //
    std::unordered_map<std::string, TemplatePosition> get_insert_size_positions(int &, std::vector<std::string>);

    // get_insert_size_positions
    //
    std::unordered_map<std::string, TemplatePosition> get_insert_size_positions(int &);

    // get_read_pairs
    //
    std::unordered_map<std::string, std::vector<BamRecord>> get_read_pairs(std::vector<std::string>);

    // get_read_pairs
    //
    std::unordered_map<std::string, std::vector<BamRecord>> get_read_pairs();

    // region_pileup
    //
    std::vector<BamRecord> region_pileup(std::vector<std::string>);

    // get_rID
    //
    int get_rID(std::string);

    // get_rName
    //
    std::string get_rName(int);

    // get_template_pos
    //
    TemplatePosition getTemplatePos(std::vector<BamRecord> &);

    // sweep_table
    //
    void sweep_table(std::unordered_map<std::string, std::vector<BamRecord>> &, BamRecord &);

    // sweep_table
    //
    void sweep_table(std::unordered_map<std::string, std::vector<BamRecord>> &, std::unordered_map<std::string, TemplatePosition> &);

    // sweep_table
    //
    void sweep_table(std::unordered_map<std::string, std::vector<BamRecord>> &, std::unordered_map<std::string, TemplatePosition> &, BamRecord &);

    // get contig lengths
    //
    std::unordered_map<std::string, int> getContigLengths();

    // get contig info
    //
    ContigInfo getContigInfo();

    // getSampleName
    // 
    // extracts sample name from the first read group in the header as string
    std::string getSampleName();

    // getReadLength()
    //
    // returns length of the next read sequence in bam file
    int getReadLength();

    // getRGNumber()
    //
    // returns the number of read groups present in the header
    int getRGNumber();
};

#endif
