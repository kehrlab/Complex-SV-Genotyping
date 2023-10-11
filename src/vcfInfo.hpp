#ifndef VCFINFOHEADER
#define VCFINFOHEADER

#include <vector>
#include <string>
#include <seqan/vcf_io.h>
#include <stdexcept>

#include "allele.hpp"
#include "genotypeResult.hpp"


class VcfInfo
{
    int rID;
    int beginPos;
    std::string id;
    std::string ref;
    std::string alt;
    int qual;
    std::string filter;
    std::vector<std::string> altAlleles;

    std::string svtype;
    std::vector<std::string> events;
    std::vector<std::string> mateids;

    std::string format;
    std::vector<std::string> genotype;
    std::vector<float> genotypeQuality;

    std::vector<std::vector<int>> genotypes;
    std::vector<std::vector<float>> gtQualities;

    std::vector<std::string> associatedAlleles;
    std::string info;

    void createAltString();
    void inferGenotypeFromAlleleCounts();
    void createInfoString();
    void setFormat();

    public:
    VcfInfo();
    VcfInfo(int, int, std::string, std::string, std::string, int, std::string, std::string, std::string, std::vector<std::string>);

    void setGenotype(std::string, float);
    void setAssociatedAlleles();
    void setRefID(int);
    void setBeginPos(int);
    void setRecordID(std::string);
    void setRefString(std::string);
    void setAltString(std::string, std::string, std::string);
    void setQuality(int);
    void setFilter(std::string);
    void setSVType(std::string);
    void setAssociatedAlleles(std::vector<std::string>);

    void extractGenotype(GenotypeResult &);

    int getRefID();
    int getBeginPos();
    std::string getRecordID();
    std::string getRefString();
    std::vector<std::string> getAltStrings();
    std::vector<std::string> getEvents();
    std::vector<std::string> getMateIDs();
    std::string getFilter();
    std::string getSVType();
    std::vector<std::vector<int>> getGenotypes();
    std::vector<std::vector<float>> getGtQualities();
    
    void mergeWith(VcfInfo);
    bool isSameBreakend(VcfInfo);
    void updateMateID(std::string, std::string);

    seqan::VcfRecord getVcfRecord();
};

#endif