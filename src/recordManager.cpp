#include "recordManager.hpp"
#include "readTemplate.hpp"
#include <unordered_map>

RecordManager::RecordManager()
{
    this->maxReadLength = 100;
}

void RecordManager::setBamRecords(std::vector<BamRecord> records)
{
    this->bamRecords = records;
    createRecordHashTable();
    collectTemplates();
}

void RecordManager::setReadPairs(std::unordered_map<std::string, std::vector<BamRecord>> recordTable)
{
    resetRecords();
    this->recordTable = recordTable;
    for (auto it : this->recordTable)
	    for (auto r : it.second)
		    this->maxReadLength = std::max(this->maxReadLength, r.getSeqLength());
    collectTemplates();
}

void RecordManager::createRecordHashTable()
{
    resetHashTable();
    for (auto it = this->bamRecords.begin(); it != this->bamRecords.end(); ++it)
    {
        this->maxReadLength = std::max(this->maxReadLength, it->getSeqLength());
        std::string hKey = std::string(it->getTemplateName());
        auto tIt = this->recordTable.find(hKey);
        if (tIt != this->recordTable.end())
        {
            this->recordTable[hKey].push_back(*it);
        } else {
            std::vector<BamRecord> tempVector;
            tempVector.push_back(*it);
            this->recordTable[hKey] = tempVector; 
        }
    }
    resetRecords();
}

void RecordManager::resetHashTable()
{
    this->recordTable.erase(this->recordTable.begin(), this->recordTable.end());
    this->recordTable = std::unordered_map<std::string, std::vector<BamRecord>>();
}

void RecordManager::collectTemplates()
{
    resetTemplates();
    for (auto it = this->recordTable.begin(); it != this->recordTable.end(); ++it)
    {
        if (it->second.size() > 1)
        {
            ReadTemplate rTemplate(it->second);
            if (rTemplate.isProperPair()) 
                this->templates.push_back(rTemplate);
        }
    }
    resetHashTable();
}

void RecordManager::resetTemplates()
{
    this->templates.erase(this->templates.begin(), this->templates.end());
    this->templates = std::vector<ReadTemplate>();
}

void RecordManager::extractInsertSizes()
{
    resetInsertSizes();
    for (auto it : this->templates)
    {
        int insertSize = it.getInsertSize();
        if (insertSize >= 0 && insertSize <= 1000) {
            this->insertSizes.push_back(insertSize);
        }
    }
}

void RecordManager::resetInsertSizes()
{
    this->insertSizes.erase(this->insertSizes.begin(), this->insertSizes.end());
    this->insertSizes = std::vector<int>();
}

std::vector<int> & RecordManager::getInsertSizes()
{
    extractInsertSizes();
    return this->insertSizes;
}

std::unordered_map<std::string, std::vector<BamRecord>> & 
RecordManager::getRecordTable()
{
    return this->recordTable;
}

std::vector<ReadTemplate> & RecordManager::getTemplates()
{
    return this->templates;
}

void RecordManager::resetRecords()
{
    this->bamRecords.erase(this->bamRecords.begin(), this->bamRecords.end());
    this->bamRecords = std::vector<BamRecord>();
}

void RecordManager::reset()
{
    resetTemplates();
    resetHashTable();
    resetInsertSizes();
    resetRecords();
}

int RecordManager::getMaxReadLength()
{
    return this->maxReadLength;
}
