#ifndef RECORDMANAGERHEADER
#define RECORDMANAGERHEADER

#include <unordered_map>
#include <vector>
#include <seqan/bam_io.h>

#include "custom_types.hpp"
#include "readTemplate.hpp"
#include "bamFileHandler.hpp"
#include "record.hpp"

class RecordManager
{
    std::vector<BamRecord> bamRecords;
    std::unordered_map<std::string, std::vector<BamRecord>> recordTable;
    std::vector<ReadTemplate> templates;

    std::vector<int> insertSizes;
    int maxReadLength;

    ContigInfo cInfo;


    void extractInsertSizes();

    void createRecordHashTable();
    void collectTemplates();

    void resetRecords();
    void resetHashTable();
    void resetTemplates();
    void resetInsertSizes();

    public:
    RecordManager(ContigInfo);

    void setBamRecords(std::vector<BamRecord>);
    void setReadPairs(std::unordered_map<std::string, std::vector<BamRecord>>);
    void reset();

    std::vector<int> & getInsertSizes();
    std::unordered_map<std::string, std::vector<BamRecord>> & getRecordTable();
    std::vector<ReadTemplate> & getTemplates();
    int getMaxReadLength();
};

#endif