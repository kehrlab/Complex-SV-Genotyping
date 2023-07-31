#ifndef RECORDMANAGERHEADER
#define RECORDMANAGERHEADER

#include <vector>
#include <seqan/bam_io.h>
#include <unordered_map>
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


    void extractInsertSizes();

    void createRecordHashTable();
    void collectTemplates();

    void resetRecords();
    void resetHashTable();
    void resetTemplates();
    void resetInsertSizes();

    public:
    RecordManager();

    void setBamRecords(std::vector<BamRecord>);
    void setReadPairs(std::unordered_map<std::string, std::vector<BamRecord>>);
    void reset();

    std::vector<int> & getInsertSizes();
    std::unordered_map<std::string, std::vector<BamRecord>> & getRecordTable();
    std::vector<ReadTemplate> & getTemplates();
    int getMaxReadLength();
};

#endif