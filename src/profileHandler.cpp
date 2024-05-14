#include "profileHandler.hpp"
#include "readTemplate.hpp"
#include "record.hpp"
#include <string>

PopDelProfileHandler::PopDelProfileHandler()
{
    this->profilePath = "";
    std::string sampleName = "";

    this->sampleNumRegions = 0;
    this->indexRegionSize = 0;
    this->numReadGroups = 0;
}


PopDelProfileHandler::PopDelProfileHandler(std::string profilePath)
{
    this->profilePath = profilePath;
    open(this->profilePath);
}


void PopDelProfileHandler::close()
{
    // close the stream
    if (this->infile.is_open())
        this->infile.close();
}


void PopDelProfileHandler::open(std::string profilePath)
{
    this->profilePath = profilePath;
    
    // open the profile
    this->infile = tryOpenHistogram(seqan::CharString(this->profilePath.c_str()));

    // load the histograms
    readHistograms();

    // merge the histograms and create sample distribution
    mergeRGHistograms();

    // extract any other relevant information, if any
}


inline void PopDelProfileHandler::readHistograms()
{
    uint16_t profileVersion;
    
    seqan::CharString sampleName;
    seqan::String<seqan::CharString> sampleReadGroups;
    seqan::String<Histogram> sampleHistograms;
    
    bool dropContigs = false; // I'm not sure about that, might need to adjust later
    
    readProfileHeader(
        this->infile,
        profileVersion,
        this->profilePath,
        sampleName,
        sampleReadGroups,
        sampleHistograms,
        this->sampleContigNames,
        this->sampleContigLengths,
        this->sampleNumRegions,
        this->indexRegionSize,
        dropContigs
    );

    // get sample name
    this->sampleName = std::string(seqan::toCString(sampleName));

    // get Histrogram for every read group
    for (unsigned i = 0; i < seqan::length(sampleReadGroups); ++i)
    {
        this->readGroups.push_back(std::string(seqan::toCString(sampleReadGroups[i])));
        this->histograms.push_back(sampleHistograms[i]);
    }
    this->numReadGroups = this->readGroups.size();
}

inline void PopDelProfileHandler::mergeRGHistograms()
{
    unsigned readLength = 0;
    int sMin = 1;
    int sMax = 1000;

    std::vector<unsigned> readLengths;

    for (unsigned i = 0; i < this->histograms.size(); ++i)
    {
        readLengths.push_back(this->histograms[i].readLength);

        // determine boundaries
        uint16_t offset = this->histograms[i].offset;
        sMin = std::min(offset + 1, sMin);
        sMax = std::max((int) (offset + seqan::length(this->histograms[i].values) - 2), sMax);
    }
    std::sort(readLengths.begin(), readLengths.end());
    readLength = readLengths[(int) (readLengths.size() / 2)];

    // add values
    std::vector<double> aggHist(sMax - sMin + 1);
    for (unsigned i = 0; i < this->histograms.size(); ++i)
    {
        uint16_t offset = this->histograms[i].offset;
        for (unsigned j = 1; j < seqan::length(this->histograms[i].values) - 1; ++j)
            aggHist[offset + j - sMin] += this->histograms[i].values[j];
    }

    // initialize library distribution from this
    this->sampleDistribution = LibraryDistribution(sMin, sMax, readLength, aggHist);

    // empty the individual histograms
    std::vector<Histogram>().swap(this->histograms);
    std::vector<std::string>().swap(this->readGroups);
}

std::unordered_map<std::string, std::vector<BamRecord>> PopDelProfileHandler::get_read_pairs(std::vector<GenomicRegion> regions)
{
    std::unordered_map<std::string, std::vector<BamRecord>> records;

    uint32_t templateIdx = 0;

    for (GenomicRegion & r : regions)
    {
        uint32_t beginPos = r.getRegionStart();
        seqan::CharString rName(r.getReferenceName().c_str());

        // jump to first window overlapping the region of interest
        jumpToRegion(this->infile, this->sampleContigNames, this->sampleContigLengths, this->indexRegionSize, rName, beginPos);
        zlib_stream::zip_istream unzipper(this->infile);

        // read alignments from windows
        Window w;
        while (readWindow(w, unzipper, this->numReadGroups))
        {
            // exit if landed on wrong chromosome
            if (this->sampleContigNames[w.chrom] != rName)
                break;

            // abort when no more windows overlap r
            if (w.beginPos > r.getRegionEnd())
                break;

            // continue if not yet at correct position
            if (w.beginPos + 255 < beginPos)
                continue;
            // continue if window is empty
            if (seqan::length(w.records) == 0)
                continue;
            
            // re-create alignment records from information and add to records
            for (unsigned rg = 0; rg < seqan::length(w.records); ++rg)
            {
                for (unsigned i = 0; i < seqan::length(w.records[rg]); ++i)
                {
                    std::vector<BamRecord> tempRecords;
                    createReadPair(tempRecords, w.records[rg][i], rName, this->sampleDistribution.getReadLength());
                    
                    if (!tempRecords[0].passesStandardFilter() || !tempRecords[1].passesStandardFilter())
                        continue;
                    
                    records[std::to_string(templateIdx)] = tempRecords;
                    ++templateIdx;
                }
            }
        }


        // jump to translocation block 
        unsigned translocIndexBeginPos = translocationIndexBeginPos(this->sampleNumRegions);
        jumpToTranslocBlock(this->infile, translocIndexBeginPos);

        
        zlib_stream::zip_istream transloc_unzipper(this->infile);
        TranslocationWindow window;

        do // no indexing in this block; read windows until the region of interest is reached and extract records
        {
            unsigned ret = readTillRoi(window, transloc_unzipper, rName, beginPos, this->sampleContigNames, this->numReadGroups);
            if (ret > 0)
                break;

            if (seqan::length(window.records) == 0)
                continue;

            // get all records from window
            for (unsigned rg = 0; rg < seqan::length(window.records); ++rg)
            {
                for (unsigned i = 0; i < seqan::length(window.records[rg]); ++i)
                {
                    std::vector<BamRecord> tempRecords;
                    createReadPair(tempRecords, window.records[rg][i], rName, this->sampleDistribution.getReadLength());

                    if (!tempRecords[0].passesStandardFilter() || !tempRecords[1].passesStandardFilter())
                        continue;
                    
                    // translocation pairs are duplicated in popdel profiles; keep only one (unless there is an accidental duplicate)
                    std::string key = tempRecords[0].getTemplateName();
                    int idx = 0;
                    bool equal {false};

                    while (records.find(key + std::to_string(idx)) != records.end())
                    {
                        equal = (records[key + std::to_string(idx)][0] == tempRecords[0] && records[key + std::to_string(idx)][1] == tempRecords[1]) || 
                                (records[key + std::to_string(idx)][0] == tempRecords[1] && records[key + std::to_string(idx)][1] == tempRecords[0]);
                        
                        if (equal)
                            break;
                        
                        ++idx;
                    }
                    if (!equal) 
                        records[key + std::to_string(idx)] = tempRecords;
                }
            }
        } while (this->sampleContigNames[window.chrom] == rName && (window.beginPos + 255) < r.getRegionEnd());
    }
    
    return records;
}


inline void PopDelProfileHandler::createReadPair(std::vector<BamRecord> & records, const ReadPair & rp, seqan::CharString rName, int readLength)
{
    std::string chromosome = seqan::toCString(rName);
    int32_t threePrimeLeft = rp.pos; 
    int32_t threePrimeRight = threePrimeLeft + rp.distance;
    
    std::vector<uint8_t> clipping = rp.clipping;

    // Reconstruct individual records from PopDel information
    // Use MapQ = 100 since no MapQ is known
    std::string orientation;
    if (rp.orientation == Orientation::FR) 
    {
        // Forward record
        records.push_back(BamRecord(chromosome, "PopDelPair", threePrimeLeft - (readLength - (clipping[0] + clipping[1])) + 1, threePrimeLeft, 100, readLength, clipping[1], clipping[0], false, true, true, false));
        // Reverse record
        records.push_back(BamRecord(chromosome, "PopDelPair", threePrimeRight, threePrimeRight + (readLength - (clipping[2] + clipping[3])) - 1, 100, readLength, clipping[3], clipping[2], true, true, false, true));
    } 
    else if (rp.orientation == Orientation::RR)
    {
        // two reverse records
        records.push_back(BamRecord(chromosome, "PopDelPair", threePrimeLeft, threePrimeLeft + (readLength - (clipping[0] + clipping[1])) - 1, 100, readLength, clipping[0], clipping[1], true, true, true, false));
        records.push_back(BamRecord(chromosome, "PopDelPair", threePrimeRight, threePrimeRight + (readLength - (clipping[2] + clipping[3])) - 1, 100, readLength, clipping[3], clipping[2], true, true, false, true));
    }
    else if (rp.orientation == Orientation::FF) 
    {
        // two forward records
        records.push_back(BamRecord(chromosome, "PopDelPair", threePrimeLeft - (readLength - (clipping[0] + clipping[1])) + 1, threePrimeLeft, 100, readLength, clipping[1], clipping[0], false, true, true, false));
        records.push_back(BamRecord(chromosome, "PopDelPair", threePrimeRight - (readLength - (clipping[2] + clipping[3])) + 1, threePrimeRight, 100, readLength, clipping[2], clipping[3], false, true, false, true));
    }
    else if (rp.orientation == Orientation::RF)
    {
        // Reverse record
        records.push_back(BamRecord(chromosome, "PopDelPair", threePrimeLeft, threePrimeLeft + (readLength - (clipping[0] + clipping[1])) - 1, 100, readLength, clipping[0], clipping[1], true, true, true, false));
        // Forward record
        records.push_back(BamRecord(chromosome, "PopDelPair", threePrimeRight - (readLength - (clipping[2] + clipping[3])) + 1, threePrimeRight, 100, readLength, clipping[2], clipping[3], false, true, false, true));
    }
    return;
}

inline void PopDelProfileHandler::createReadPair(std::vector<BamRecord> & records, const TranslocationWindowEntry & windowEntry, seqan::CharString rName, int readLength)
{
    std::string chromosome1 = seqan::toCString(rName);
    int32_t pos1 = windowEntry.first.pos;
    std::string chromosome2 = seqan::toCString(this->sampleContigNames[windowEntry.second.refID]);
    int32_t pos2 = windowEntry.second.pos;

    
    std::string templateName = std::min(chromosome2, chromosome1) + std::to_string(std::min(pos1, pos2)) + std::max(chromosome2, chromosome1) + std::to_string(std::max(pos1, pos2));

    if (windowEntry.orientation == Orientation::RF || windowEntry.orientation == Orientation::RR)        // first record is in reverse orientation
        records.push_back(BamRecord(chromosome1, templateName, pos1, pos1 + readLength - (windowEntry.first.clip_0 + windowEntry.first.clip_1), 100, readLength, windowEntry.first.clip_0, windowEntry.first.clip_1, true, true, true, false));
    else                                                                                                 // first record is in forward orientation
        records.push_back(BamRecord(chromosome1, templateName, pos1 - readLength + (windowEntry.first.clip_0 + windowEntry.first.clip_1), pos1, 100, readLength, windowEntry.first.clip_1, windowEntry.first.clip_0, false, true, true, false));

    if (windowEntry.orientation == Orientation::FR || windowEntry.orientation == Orientation::RR)        // second record is in reverse orientation
        records.push_back(BamRecord(chromosome2, templateName, pos2, pos2 + readLength - (windowEntry.second.clip_0 + windowEntry.second.clip_1), 100, readLength, windowEntry.second.clip_0, windowEntry.second.clip_1, true, true, false, true));
    else                                                                                                 // second record is in forward orientation
        records.push_back(BamRecord(chromosome2, templateName, pos2 - readLength + (windowEntry.second.clip_0 + windowEntry.second.clip_1), pos2, 100, readLength, windowEntry.second.clip_1, windowEntry.second.clip_0, false, true, false, true));
}


LibraryDistribution & PopDelProfileHandler::getLibraryDistribution()
{
    return this->sampleDistribution;
}

std::string PopDelProfileHandler::getSampleName()
{
    return this->sampleName;
}

std::vector<int> PopDelProfileHandler::getContigLengths()
{
    std::vector<int> contigLengths;

    for (unsigned i = 0; i < seqan::length(this->sampleContigLengths); ++i)
        contigLengths.push_back(static_cast<int>(this->sampleContigLengths[i]));

    return contigLengths;
}

std::vector<std::string> PopDelProfileHandler::getContigNames()
{
    std::vector<std::string> contigNames;

    for (unsigned i = 0; i < seqan::length(this->sampleContigNames); ++i)
        contigNames.push_back(std::string(seqan::toCString(this->sampleContigNames[i])));

    return contigNames;
}

std::string PopDelProfileHandler::getFileName()
{
    return this->profilePath;
}