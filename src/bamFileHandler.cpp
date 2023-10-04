#include "bamFileHandler.hpp"
#include "custom_types.hpp"
#include "genomicRegion.hpp"
#include "readTemplate.hpp"
#include <stdexcept>
#include <unordered_map>

BamFileHandler::BamFileHandler()
{
    this->bamFileName = "";
    this->rID = 0;
    this->beginPos = 0;
    this->endPos = 0;
    this->minQ = 0;
    this->hasAlignments = false;
}


BamFileHandler::BamFileHandler(std::string bamFileName)
{
    this->bamFileName = bamFileName;
    this->rID = 0;
    this->beginPos = 0;
    this->endPos = 0;
    this->hasAlignments = false;
    this->minQ = 0;
    inferBaiFileName();
    openInputFiles();
    extractContigInfo();
}


BamFileHandler::BamFileHandler(std::string bamFileName, ProgramOptions & options)
{
    this->bamFileName = bamFileName;
    this->rID = 0;
    this->beginPos = 0;
    this->endPos = 0;
    this->hasAlignments = false;
    this->minQ = options.getMinimumMappingQuality();
    inferBaiFileName();
    openInputFiles();
    extractContigInfo();
}


void BamFileHandler::open(std::string bamFileName, ProgramOptions & options)
{
    this->bamFileName = bamFileName;
    this->minQ = options.getMinimumMappingQuality();
    inferBaiFileName();
    openInputFiles();
    extractContigInfo();
}


void BamFileHandler::inferBaiFileName()
{
    this->baiFileName = this->bamFileName;
    this->baiFileName.append(".bai");
}


void BamFileHandler::openInputFiles()
{
    openBamFile();
    openBaiFile();
    readHeader();
}


void BamFileHandler::openBamFile()
{
    bamFileIn = hts_open(this->bamFileName.c_str(), "r");
}


void BamFileHandler::openBaiFile()
{
    index = sam_index_load(this->bamFileIn, this->baiFileName.c_str());
}


void BamFileHandler::readHeader()
{
    this->header = sam_hdr_read(this->bamFileIn);
}


void BamFileHandler::closeInputFile()
{
    bam_hdr_destroy(this->header);
    hts_idx_destroy(this->index);
    sam_close(this->bamFileIn);
    resetRegions();
}


std::string BamFileHandler::getFileName()
{
    return this->bamFileName;
}


void BamFileHandler::setRegions(std::vector<GenomicRegion> regions)
{
    resetRegions();
    for (GenomicRegion & r : regions)
        this->regions.push_back(r.getRegionString());
}


void BamFileHandler::resetRegions()
{
    std::vector<std::string>().swap(this->regions);
}


std::unordered_map<std::string, TemplatePosition> BamFileHandler::get_insert_size_positions(int & maxReadLength, std::vector<std::string> regions)
{
    std::unordered_map<std::string, std::vector<BamRecord>> recordTable;
    std::unordered_map<std::string, TemplatePosition> positions;
    bam1_t * record = bam_init1();

    int readCounter = 0;
    maxReadLength = 0;
    for (auto & r : regions)
    {
        hts_itr_t * itr = sam_itr_querys(this->index, this->header, r.c_str());

        while (sam_itr_next(this->bamFileIn, itr, record) >= 0) 
        {
            if (record->core.tid < 0)
                break;
	    
            BamRecord bamRcrd(record, this->header);
            maxReadLength = std::max(maxReadLength, bamRcrd.getSeqLength());
            if (bamRcrd.passesInsertFilter())
            {
                std::string hKey = bamRcrd.getTemplateName();
                if (recordTable.find(hKey) == recordTable.end())
                    recordTable[hKey] = std::vector<BamRecord>();
                recordTable[hKey].push_back(bamRcrd);
                ++readCounter;

                if (readCounter == 10000) {
                    sweep_table(recordTable, positions, bamRcrd);
                    readCounter = 0;
                }
            }
        }
        sweep_table(recordTable, positions);
        recordTable = std::unordered_map<std::string, std::vector<BamRecord>>();
        hts_itr_destroy(itr);
    }
    bam_destroy1(record);
    return positions;
}


std::unordered_map<std::string, TemplatePosition> BamFileHandler::get_insert_size_positions(int & maxReadLength)
{
    std::unordered_map<std::string, std::vector<BamRecord>> recordTable;
    std::unordered_map<std::string, TemplatePosition> positions;
    bam1_t * record = bam_init1();
    hts_itr_t * itr = sam_itr_querys(this->index, this->header, ".");
    maxReadLength = 0;

    int readCounter = 0;
    while (sam_itr_next(this->bamFileIn, itr, record) >= 0) 
    {
        if (record->core.tid < 0)
            break;
        BamRecord bamRcrd(record, this->header);
        maxReadLength = std::max(maxReadLength, bamRcrd.getSeqLength());
        if (bamRcrd.passesInsertFilter())
        {
            std::string hKey = bamRcrd.getTemplateName();
            if (recordTable.find(hKey) == recordTable.end())
                recordTable[hKey] = std::vector<BamRecord>();
            recordTable[hKey].push_back(bamRcrd);
            ++readCounter;

            if (readCounter == 10000) {
                sweep_table(recordTable, positions, bamRcrd);
                readCounter = 0;
            }
        }
    }

    sweep_table(recordTable, positions);
    recordTable = std::unordered_map<std::string, std::vector<BamRecord>>();
    hts_itr_destroy(itr);
    bam_destroy1(record);

    return positions;
}

std::unordered_map<std::string, std::vector<BamRecord>> BamFileHandler::get_read_pairs(std::vector<std::string> regions)
{
    std::unordered_map<std::string, std::vector<BamRecord>> recordTable;
    bam1_t * record = bam_init1();

    
    for (auto & r : regions)
    {
        hts_itr_t * itr = sam_itr_querys(this->index, this->header, r.c_str());

	
        while (sam_itr_next(this->bamFileIn, itr, record) >= 0) 
        {
            if (record->core.tid < 0)
                break;
	
            BamRecord bamRcrd(record, this->header);
            if (bamRcrd.passesStandardFilter())
            {
                std::string hKey = bamRcrd.getTemplateName();
                if (recordTable.find(hKey) == recordTable.end())
                    recordTable[hKey] = std::vector<BamRecord>();
                recordTable[hKey].push_back(bamRcrd);
            }

        }
        hts_itr_destroy(itr);
    }
    bam_destroy1(record);
    
    return recordTable;
}


std::unordered_map<std::string, std::vector<BamRecord>> BamFileHandler::get_read_pairs()
{
    std::unordered_map<std::string, std::vector<BamRecord>> recordTable;
    bam1_t * record = bam_init1();

    for (auto & r : this->regions)
    {
        hts_itr_t * itr = sam_itr_querys(this->index, this->header, r.c_str());
        while (sam_itr_next(this->bamFileIn, itr, record) >= 0) 
        {
            if (record->core.tid < 0)
                break;
            BamRecord bamRcrd(record, this->header);
            if (bamRcrd.passesStandardFilter())
            {
                std::string hKey = bamRcrd.getTemplateName();
                if (recordTable.find(hKey) == recordTable.end())
                    recordTable[hKey] = std::vector<BamRecord>();
                recordTable[hKey].push_back(bamRcrd);
            }
        }
        hts_itr_destroy(itr);
    }
    bam_destroy1(record);
    
    return recordTable;
}


std::vector<BamRecord> BamFileHandler::region_pileup(std::vector<std::string> regions)
{
    std::vector<BamRecord> records;
    bam1_t * record = bam_init1();

    for (auto & r : regions)
    {
        hts_itr_t * itr = sam_itr_querys(this->index, this->header, r.c_str());
        
        while (sam_itr_next(this->bamFileIn, itr, record) >= 0)
        {
            if (record->core.tid < 0)
                break;
            BamRecord bamRcrd(record, this->header);
            if (bamRcrd.passesStandardFilter())
                records.push_back(bamRcrd);
        }
        hts_itr_destroy(itr);
    }
    bam_destroy1(record);
    return records;
}


TemplatePosition BamFileHandler::getTemplatePos(std::vector<BamRecord> & records)
{
    TemplatePosition pos;
    pos.begin = 0;
    pos.end = 0;
    pos.chr = "";

    if (records[0].getReferenceName() != records[1].getReferenceName())
        return pos;
    pos.chr = records[0].getReferenceName();
    pos.begin = std::min(records[0].getFivePrimePos(), records[1].getFivePrimePos());
    pos.end = std::max(records[0].getFivePrimePos(), records[1].getFivePrimePos());
    return pos;
}

void BamFileHandler::sweep_table(std::unordered_map<std::string, std::vector<BamRecord>> & recordTable, std::unordered_map<std::string, TemplatePosition> & positionTable)
{
    for (auto it = recordTable.begin(); it != recordTable.end(); ++it)
    {
        if (it->second.size() == 2)
        {
            int leftIdx, rightIdx;
            if (it->second[0].getStartPos() < it->second[1].getStartPos())
            {
                leftIdx = 0;
                rightIdx = 1;
            } else {
                leftIdx = 1;
                rightIdx = 0;
            }
            if (!it->second[leftIdx].isReverse() && it->second[rightIdx].isReverse())
            {
                TemplatePosition tempPos = getTemplatePos(it->second);
                int is = tempPos.end - tempPos.begin + 1;
                if (is < 1000 && is > 0 && tempPos.chr != "")
                    positionTable[it->first] = tempPos;
            }
        }
    }
}


void BamFileHandler::sweep_table(std::unordered_map<std::string, std::vector<BamRecord>> & recordTable, std::unordered_map<std::string, TemplatePosition> & positionTable, BamRecord & bamRcrd)
{
    for (auto it = recordTable.begin(); it != recordTable.end(); )
    {
        if (it->second.size() == 2)
        {
            int leftIdx, rightIdx;
            if (it->second[0].getStartPos() < it->second[1].getStartPos())
            {
                leftIdx = 0;
                rightIdx = 1;
            } else {
                leftIdx = 1;
                rightIdx = 0;
            }
            if (!it->second[leftIdx].isReverse() && it->second[rightIdx].isReverse())
            {
                TemplatePosition tempPos = getTemplatePos(it->second);
                int is = tempPos.end - tempPos.begin + 1;
                if (is < 1000 && is > 0 && tempPos.chr != "")
                    positionTable[it->first] = tempPos;
            }
            it = recordTable.erase(it);
            continue;
        }
        if (it->second.size() > 2)
        {
            it = recordTable.erase(it);
            continue;
        }
        if ((it->second[0].getStartPos() - bamRcrd.getEndPos()) > 1000) {
            it = recordTable.erase(it);
            continue;
        }
        ++it;
    }
}


int BamFileHandler::get_rID(std::string refName)
{
    return sam_hdr_name2tid(this->header, refName.c_str());
}


std::string BamFileHandler::get_rName(int rID)
{
    return sam_hdr_tid2name(this->header, rID);
}


std::unordered_map<std::string, int> BamFileHandler::getContigLengths()
{
    std::unordered_map<std::string, int> contigLengths;

    for (int i = 0; i < this->contigInfo.cNames.size(); ++i)
    {
        std::string rName = this->contigInfo.cNames[i];
        int64_t rLen = this->contigInfo.cLengths[i];
        contigLengths[rName] = rLen;
    }
    return contigLengths;
}


void BamFileHandler::extractContigInfo()
{
    ContigInfo cInfo;
    int nref = sam_hdr_nref(this->header);
    for (int i = 0; i < nref; ++i)
    {
        cInfo.cNames.push_back(get_rName(i));
        cInfo.cLengths.push_back(sam_hdr_tid2len(this->header, i));
    }
    this->contigInfo = cInfo;
}


ContigInfo BamFileHandler::getContigInfo()
{
    return this->contigInfo;
}
