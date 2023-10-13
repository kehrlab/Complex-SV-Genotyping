#ifndef VCFWRITERHEADER
#define VCFWRITERHEADER

#include <algorithm>
#include <ctime>
#include <stdexcept>
#include <string>
#include <vector>
#include <seqan/vcf_io.h>
#include <ctime>

#include "bamFileHandler.hpp"
#include "variant.hpp"
#include "genotypeResult.hpp"
#include "vcfInfo.hpp"


class VcfWriter
{
    std::string vcfFileName;
    seqan::VcfFileOut vcfFile;
    seqan::VcfHeader header;

    std::vector<std::string> sampleNames;
    std::unordered_set<std::string> contigNames;
    std::unordered_map<std::string, int> rIDs;

    std::vector<seqan::VcfRecord> vcfRecords;
    
    void writeHeader();
    void setFormalities();
    void setSampleNames();
    void setContigNames();
    void initHeader();
    void createRefIDs();

    void writeVcfRecords();
    void joinRecords(seqan::VcfRecord &, seqan::VcfRecord);

    public:
    VcfWriter();
    VcfWriter(std::string);

    void addVariantRecords(std::vector<GenotypeResult> &, complexVariant &);
    void setFileName(std::string);
    void openVcfFile();
    void closeVcfFile();
    void write();
};

#endif