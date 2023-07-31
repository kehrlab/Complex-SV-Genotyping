#ifndef VCFWRITERHEADER
#define VCFWRITERHEADER

#include "bamFileHandler.hpp"
#include "variant.hpp"
#include <string>
#include <vector>
#include <seqan/vcf_io.h>
#include <ctime>
#include "genotypeResult.hpp"

class VcfWriter
{
    std::string vcfFileName;
    seqan::VcfFileOut vcfFile;
    seqan::VcfHeader header;

    std::vector<std::string> sampleNames;
    std::vector<std::string> contigNames;
    std::unordered_map<std::string, int> rIDs;

    std::vector<seqan::VcfRecord> vcfRecords;
    
    void getContigNames(std::vector<complexVariant> &);
    void getSampleNames(std::vector<std::string> fileNames);
    void writeHeader();
    void setFormalities();
    void setSampleNames();
    void setContigNames();
    void initHeader();
    void createRefIDs();

    void createVcfRecords(std::vector<complexVariant> &, std::vector<std::vector<GenotypeResult>> &);
    void writeVcfRecords();
    void joinRecords(seqan::VcfRecord &, seqan::VcfRecord);

    public:
    VcfWriter();
    VcfWriter(std::string);
    VcfWriter(std::string, std::vector<complexVariant> &, std::vector<std::string>, std::vector<std::vector<GenotypeResult>> &);

    void setFileName(std::string);
    void openVcfFile();
    void closeVcfFile();
    void write();
};

#endif