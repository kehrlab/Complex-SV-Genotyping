#ifndef POPDELSNIPPETSH
#define POPDELSNIPPETSH

#include <cstdint>
#include <seqan/sequence.h>
#include <seqan/file/file_interface.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>

// Code snippets taken from PopDel with only minor changes

// -----------------------------------------------------------------------------
// Histogram
// -----------------------------------------------------------------------------


// -----------------------------------------------------------------------------
// Typedefs
// -----------------------------------------------------------------------------

typedef seqan::Iterator<seqan::String<double>, seqan::Rooted>::Type THistIter;
// -----------------------------------------------------------------------------
// Struct Histogram
// -----------------------------------------------------------------------------
struct Histogram
{                          // Index 1 holds the counts of distances of len. 1 (+offset) etc
    seqan::String<double> values; // First and last index are reserved for irregular values (too big, etc.)
    double         min_prob;               // 1/10000 of the maximum count in "values"
    uint16_t       offset;                 // Difference of the first element in the string from inner distance 0
    double         mean;
    double         stddev;
    uint16_t       median;                 // The median insert size
    uint16_t       median3PrimeDist;       // The median distance between the 3'ends. Values below 0 are set to 0
    unsigned       readLength;             // length of one read making up the read pair
    double         coverage;
    unsigned       lowerQuantileDist;      // Distance from median to 1% quantile.
    unsigned       upperQuantileDist;      // Distance from median to 99% quantile.
    unsigned       windowSize;             // The size of one window in bp's.

    Histogram():
    values(""),
    min_prob(0.0),
    offset(0),
    mean(0.0),
    stddev(0.0),
    median(0),
    median3PrimeDist(0),
    readLength(0),
    coverage(0),
    windowSize(0)
    {}
    // ---------------------------------------------------------------------------------------
    // Function clearStrings()
    // ---------------------------------------------------------------------------------------
    // Clears all member strings. If shrink is true, the object will be shrunk to capacity 0.
    inline void clearStrings(bool shrink=false)
    {
        clear(values);
        if (shrink)
        {
            resize(values, 0);
            shrinkToFit(values);
        }
        else
        {
            clear(values);
        }
    }
    // Return the size of the histogram in byte
    inline unsigned getSize() const
    {
        return (sizeof(values) +
                length(values) * sizeof(double) +
                4 * sizeof(double) +
                1 * sizeof(int) +
                5 * sizeof(unsigned));
    }
    // Return the capcity of the histogram in byte
    inline unsigned getCapactiy() const
    {
        return (sizeof(values) +
                capacity(values) * sizeof(double) +
                4 * sizeof(double) +
                1 * sizeof(int) +
                5 * sizeof(unsigned));
    }
};


// -----------------------------------------------------------------------------
// Function tryOpenHistogram()
// -----------------------------------------------------------------------------
// Trie to open the histogram file and throw an error on failures.
// Return the opened ifstream.
inline std::ifstream tryOpenHistogram(const seqan::CharString & filename)
{
    std::ifstream infile(toCString(filename));
    if (!infile.is_open())
    {
        std::ostringstream msg;
        msg << "[PopDel] Could not open histogram file \'" << filename << "\' for reading.";
        SEQAN_THROW(seqan::IOError(seqan::toCString(msg.str())));
    }
    return infile;
}

// =======================================================================================
// Function readProfileHeader()
// =======================================================================================
template<typename TStream>
inline void readProfileHeader(TStream & stream,
                              uint16_t & version,
                              const seqan::CharString & filename,
                              seqan::CharString & sampleName,
                              seqan::String<seqan::CharString> & readGroups,
                              seqan::String<Histogram> & histograms,
                              seqan::String<seqan::CharString> & contigNames,
                              seqan::String<int32_t> & contigLengths,
                              unsigned & numRegions,
                              unsigned & indexRegionSize,
                              bool dropContigs = false)
{
    // Read the magic string.
    seqan::CharString buffer;
    resize(buffer, 7);
    stream.read(&buffer[0], 7);

    //Prepare message for potential errors.
    std::ostringstream msg;
    msg << "[PopDel] Corrupted profile \"" << filename << "\": ";

    if (!stream.good())
    {
        msg << "Unable to read magic string.";
        SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
    }
    if (buffer != "POPDEL\1")
    {
        msg << "Magic string is wrong.";
        SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
    }
    // Read the version
    stream.read(reinterpret_cast<char *>(&version), sizeof(uint16_t));
    if (!stream.good())
    {
        msg << "Unable to read profile version.";
        SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
    }
    // Read the index regions size, size of the index and skip (ignore) the index.
    stream.read(reinterpret_cast<char *>(&indexRegionSize), sizeof(uint32_t));
    if (!stream.good())
    {
        msg << "Unable to read index region size.";
        SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
    }
    stream.read(reinterpret_cast<char *>(&numRegions), sizeof(uint32_t));
    if (!stream.good())
    {
        msg << "Unable to read index size.";
        SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
    }
    unsigned numTranslocRegions = 0;
    stream.read(reinterpret_cast<char *>(&numTranslocRegions), sizeof(uint32_t));
    if (!stream.good())
    {
        msg << "Unable to read translocation index size.";
        SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
    }
    stream.ignore((numRegions + numTranslocRegions) * sizeof(uint64_t));

    // Read the length of the sample name.
    uint32_t smLen = 0;
    stream.read(reinterpret_cast<char *>(&smLen), sizeof(uint32_t));
    if (!stream.good())
    {
        msg << "Unable to read length of sample name.";
        SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
    }
    // Read the sample name
    resize(sampleName, smLen-1, seqan::Exact());
    stream.read(reinterpret_cast<char *>(&sampleName[0]), smLen-1);
    stream.read(&buffer[0], 1);
    if (!stream.good())
    {
        msg << "Unable to read sample name.";
        SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
    }
    if (buffer[0] != '\0') // Expect sample name to be null-terminated.
    {
        msg << "Expecting sample name to be null-terminated.";
        SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
    }
    // Read the number of read groups.
    uint32_t numReadGroups = 0;
    stream.read(reinterpret_cast<char *>(&numReadGroups), sizeof(uint32_t));
    if (!stream.good())
    {
        msg << "Unable to read number of read groups.";
        SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
    }
    if (numReadGroups == 0) // How can this even happen?
    {
        msg << "Profile seems to have 0 read groups.";
        SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
    }
    resize(readGroups, numReadGroups, seqan::Exact());
    resize(histograms, numReadGroups, seqan::Exact());

    // Read the read group names and histograms.
    for (unsigned i = 0; i < numReadGroups; ++i)
    {
        // Read length of read group name.
        uint32_t nameLen = 0;
        stream.read(reinterpret_cast<char *>(&nameLen), sizeof(uint32_t));
        if (!stream.good())
        {
            msg << "Unable to read length of read group name.";
            SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
        }

        // Read the read group name.
        resize(readGroups[i], nameLen-1, seqan::Exact());
        stream.read(reinterpret_cast<char *>(&readGroups[i][0]), nameLen-1);
        stream.read(&buffer[0], 1);
        if (!stream.good())
        {
            msg << "Unable to read read group name.";
            SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
        }
        if (buffer[0] != '\0') // Expect read group names to be null-terminated.
        {
            msg << "Expecting read group name to be null-terminated.";
            SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
        }
        // Read the statistics.
        stream.read(reinterpret_cast<char *>(&histograms[i].median), sizeof(uint16_t));
        if (!stream.good())
        {
            msg << "Unable to read median.";
            SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
        }
        stream.read(reinterpret_cast<char *>(&histograms[i].stddev), sizeof(double));
        if (!stream.good())
        {
            msg << "Unable to read standard deviation of distances";
            SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
        }
        stream.read(reinterpret_cast<char *>(&histograms[i].readLength), sizeof(uint16_t));
        if (!stream.good())
        {
            msg << "Unable to read read length.";
            SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
        }
        // calculateMedian3PrimeDist(histograms[i]);
        stream.read(reinterpret_cast<char *>(&histograms[i].offset), sizeof(uint16_t));
        if (!stream.good())
        {
            msg << "Unable to read histogram offset.";
            SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
        }
        unsigned histEnd = 0;
        stream.read(reinterpret_cast<char *>(&histEnd), sizeof(uint16_t));
        if (!stream.good())
        {
            msg << "Unable to read histogram size.";
            SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
        }
        resize(histograms[i].values, histEnd - histograms[i].offset, seqan::Exact());
        for (unsigned h = 0; h < histEnd - histograms[i].offset; ++h)
        {
            stream.read(reinterpret_cast<char *>(&histograms[i].values[h]), sizeof(double));
            if (!stream.good())
            {
                msg << "Unable to read value from histogram.";
                SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
            }
        }
    }

    // Read the number of chromosomes.
    uint32_t numContigs = 0;
    stream.read(reinterpret_cast<char *>(&numContigs), sizeof(uint32_t));
    if (!stream.good())
    {
        msg << "Unable to read number of contigs.";
        SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
    }
    if (!dropContigs)
    {
        resize(contigNames, numContigs, seqan::Exact());
        resize(contigLengths, numContigs, seqan::Exact());
    }
    // Read the chromosome names and lengths.
    for (unsigned i = 0; i < numContigs; ++i)
    {
        // Read length of chromosome name.
        uint32_t nameLen = 0;
        stream.read(reinterpret_cast<char *>(&nameLen), sizeof(uint32_t));
        if (!stream.good())
        {
            msg << "Unable to read length of contig name.";
            SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
        }

        // Read the chromosome name.
        if (!dropContigs)
        {
            resize(contigNames[i], nameLen-1, seqan::Exact());
            stream.read(reinterpret_cast<char *>(&contigNames[i][0]), nameLen-1);
            stream.read(&buffer[0], 1);
            if (!stream.good())
            {
                msg << "Unable to read contig name.";
                SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
            }
            if (buffer[0] != '\0') // Expect chromosome names to be null-terminated.
            {
                msg << "Expecting contig name to be null-terminated.";
                SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
            }

            // Read the chromosome length.
            stream.read(reinterpret_cast<char *>(&contigLengths[i]), sizeof(int32_t));
            if (!stream.good())
            {
                msg << "Unable to read contig length.";
                SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
            }
        }
        else
        {
            stream.ignore(nameLen + sizeof(int32_t));
            if (!stream.good())
            {
                msg << "Unable to skip contig name and length.";
                SEQAN_THROW(seqan::ParseError(seqan::toCString(msg.str())));
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Window
// -----------------------------------------------------------------------------

// =======================================================================================
// Enum Struct Orientation
// =======================================================================================
// Enum for orientation of the read pairs.
enum struct Orientation: uint8_t
{
    FR = 0,     // > <
    FF = 1,     // > >
    RR = 2,     // < <
    RF = 3      // < >
};
// =======================================================================================
// Function firstReadForward()
// =======================================================================================
// Return true if the first read of the read pair is pointing forward.
// Return false otherwise
inline bool firstReadForward(const Orientation & o)
{
    return o == Orientation::FR || o == Orientation::FF;
}
// =======================================================================================
// Function secondReadForward()
// =======================================================================================
// Return true if the second read of the read pair is pointing forward.
// Return false otherwise
inline bool secondReadForward(const Orientation & o)
{
    return o == Orientation::RF || o == Orientation::FF;
}
// =======================================================================================
// Function getPairOrientation()
// =======================================================================================
// Return the orientation of the read pair.
// IMPORTANT: Assumes that the given record is the rightmost read of the pair!
inline Orientation getPairOrientation(const uint32_t flag)
{
    if (flag & 16)           // read reverse strand
    {
        if (flag & 32)       // mate reverse strand
        {
            return Orientation::RR;
        }
        else
        {
            return  Orientation::FR;
        }
    }
    else if (flag & 32)      // mate reverse strand
    {
        return Orientation::RF;
    }
    else
    {
        return Orientation::FF;
    }
}
// Overload for directly applying getPairOrientation() on a BamAlignmentRecord.
inline Orientation getPairOrientation(const seqan::BamAlignmentRecord & record)
{
    return getPairOrientation(record.flag);
}

// std::ostream & operator << (std::ostream & out, const Orientation & o)
// {
//     if (o == Orientation::FR )
//         out << "FR";
//     else if (o == Orientation::FF )
//         out << "FF";
//     else if (o == Orientation::RR )
//         out << "RR";
//     else
//         out << "RF";
//     return out;
// }

// =======================================================================================
// Function getSwitchedOrientation()
// =======================================================================================
// Return the orientation one would get by switching first and second read, but keeping their individual orientations.
// FR -> RF, RF -> FR, but FF -> FF and RR -> RR.
inline Orientation getSwitchedOrientation(const Orientation o)
{
    if (o == Orientation::FR)
        return Orientation::RF;
    else if (o == Orientation::RF)
        return Orientation::FR;
    else
        return o;
}
// =======================================================================================
// Struct ReadPair
// =======================================================================================
// Holds the properties of one record (=read pair) in the window.
struct ReadPair
{
    unsigned pos;                   // Position of the 3' end of the leftmost read in the pair
    int32_t distance;                   // FR-reads: 5'-end distance. Other reads; 3'-end distance
    std::vector<uint8_t> clipping;   // number of clipped bases (5', left; 3', left; 3', right; 5', right)
    uint16_t totalClipping;
    Orientation orientation;

    ReadPair():
        pos(0), distance(0), clipping(std::vector<uint8_t>(4, 0)), totalClipping(0), orientation(Orientation::FR )
    {}

    ReadPair(unsigned p, unsigned i, uint8_t c, Orientation o = Orientation::FR ):
        pos(p), distance(i), clipping({0, c/2, c/2, 0}), totalClipping(c), orientation(o)
    {}

    ReadPair(unsigned p, unsigned i, uint8_t c_0, uint8_t c_1, uint8_t c_2, uint8_t c_3, Orientation o = Orientation::FR ):
        pos(p), distance(i), clipping({c_0, c_1, c_2, c_3}), totalClipping(c_0 + c_1 + c_2 + c_3), orientation(o)
    {}

    ReadPair(unsigned p, unsigned i, std::vector<uint8_t> c, Orientation o = Orientation::FR ):
        pos(p), distance(i), clipping(c), totalClipping(0), orientation(o)
    {
        for (unsigned j = 0; j < clipping.size(); ++j)
            totalClipping += clipping[j];
    }

    bool operator <(const ReadPair & r)
    {
        if (pos < r.pos)
            return true;
        if (pos == r.pos && distance < r.distance )
            return true;
        else
            return false;
    }
};
inline bool operator==(const ReadPair& l, const ReadPair& r)
{
    return l.pos == r.pos && l.distance == r.distance && l.clipping[0] == r.clipping[0] && l.clipping[1] == r.clipping[1] && l.clipping[2] == r.clipping[2] && l.clipping[3] == r.clipping[3] && l.orientation == r.orientation;
}
inline bool operator!=(const ReadPair& l, const ReadPair& r)
{
    return !(l == r);
}
// =======================================================================================
// Struct Window
// =======================================================================================
// Defines a window on the Genome. Can contain samples from multiple read groups.
struct Window
{
    int32_t chrom;
    int32_t beginPos;
    seqan::String<seqan::String<ReadPair> > records;

    Window():
        chrom(seqan::GenomicRegion::INVALID_ID), beginPos(seqan::GenomicRegion::INVALID_POS)
    {}

    Window(__int32 c, __int32 pos):
        chrom(c), beginPos(pos)
    {}
};

// =======================================================================================
// Function reset()
// =======================================================================================
// Reset chrom, beginPos and records window but keep the same number of read groups.
template<typename TWindow>
inline void reset(TWindow & window)
{
    window.chrom = 0;
    window.beginPos = 0;
    for (unsigned rg = 0; rg < seqan::length(window.records); ++rg)
        seqan::clear(window.records[rg]);
}

// =======================================================================================
// Function readCoordinates()
// =======================================================================================
// Read the chromosome and position of the window from the input stream. Used by readWindow().
// Return true on success, false if EOF has been reached or the translocation guard is encountered
inline bool readCoordinates(int32_t & chrom, int32_t & beginPos, zlib_stream::zip_istream & stream)
{
    // Read chromosome and position.
    uint32_t refID = chrom;
    stream.read(reinterpret_cast<char *>(&refID), sizeof(uint32_t));
    if (refID == seqan::maxValue<uint32_t>() || stream.eof())
        return false;
    SEQAN_ASSERT_LEQ(refID, static_cast<uint32_t>(seqan::maxValue<int32_t>()));
    chrom = refID;
    if (!stream.good())
        SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read contig name."));
    stream.read(reinterpret_cast<char *>(&beginPos), sizeof(uint32_t));
    if (!stream.good())
        SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read contig position."));
    return true;
}

// =======================================================================================
// Function readNumRecords()
// =======================================================================================
// Return the number of read pairs in the current window of the input stream. Used by readWindow().
inline unsigned readNumRecords(zlib_stream::zip_istream & stream)
{
    // Read the number of read pairs for this read group.
    uint32_t numRecords = 0;
    stream.read(reinterpret_cast<char *>(&numRecords), sizeof(uint32_t));
    if (!stream.good())
        SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read number of read pairs."));
    return numRecords;
}

// =======================================================================================
// Function readNumClippedRecords()
// =======================================================================================
// Return the number of clipped read pairs in the current window of the input stream. Used by readWindow().
inline unsigned readNumClippedRecords(zlib_stream::zip_istream & stream)
{
    // Read the number of read pairs for this read group.
    uint32_t numClippedRecords = 0;
    stream.read(reinterpret_cast<char *>(&numClippedRecords), sizeof(uint32_t));
    if (!stream.good())
        SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read number of clipped read pairs."));
    return numClippedRecords;
}

// =======================================================================================
// Function readNumClippedRecords()
// =======================================================================================
// Return the number of clipped read pairs in the current window of the input stream. Used by readWindow().
inline unsigned readNumNonFRRecords(zlib_stream::zip_istream & stream)
{
    // Read the number of read pairs for this read group.
    uint32_t numNonFRRecords = 0;
    stream.read(reinterpret_cast<char *>(&numNonFRRecords), sizeof(uint32_t));
    if (!stream.good())
        SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read number non-FR oriented records."));
    return numNonFRRecords;
}

// =======================================================================================
// Function unpackOrientations()
// =======================================================================================
// Put the orientations in groups of 4, s.t. they can fill all 8 bits of a uint8_t. Put these uints8_t into a String.
inline void unpackOrientations(seqan::String<ReadPair> & records, const seqan::String<uint8_t> & packed, const seqan::String<unsigned> & ids)
{
    uint8_t mask = 192;     // Translates to 11000000 in binary format.
    for (unsigned i = 0; i < seqan::length(ids); ++i)
        records[ids[i]].orientation = static_cast<Orientation>((packed[i / 4] & (mask >> i % 4 * 2)) >> (6-i%4*2));
}

// =======================================================================================
// Function readOrientations()
// =======================================================================================
// Read the packed string of read pair orientations from the stream and unpack them into the records.
// Used by readWindow().
inline void readOrientations(seqan::String<ReadPair> & records, zlib_stream::zip_istream & stream, const unsigned & numRecords)
{
    if (numRecords == 0)
        return;

    seqan::String<uint8_t> orientations;
    seqan::String<uint32_t> ids;
    unsigned orientationChars = std::ceil(static_cast<double>(numRecords) / 4);
    seqan::resize(orientations, orientationChars, seqan::Generous());
    seqan::resize(ids, numRecords, seqan::Exact());
    stream.read(reinterpret_cast<char *>(&ids[0]), sizeof(uint32_t) * numRecords);
    if (!stream.good())
        SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read orientation IDs of read pairs."));
    stream.read(reinterpret_cast<char *>(&orientations[0]), orientationChars);
    if (!stream.good())
        SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read orientations of read pairs."));
    unpackOrientations(records, orientations, ids);   // Unpack the orientations and store them in the window.
}

// =======================================================================================
// Function readClipping()
// =======================================================================================
// Read all clipping info of the window's records and create a map of ID:clipping.
inline void readClipping(std::map<unsigned, std::vector<uint8_t>> & clipMap,
                         zlib_stream::zip_istream & stream,
                         const unsigned & numClippedRecords)
{
    for (unsigned i = 0; i < numClippedRecords; ++i)
    {
        unsigned id = 0;
        std::vector<uint8_t> clipping(4, 0);

        stream.read(reinterpret_cast<char *>(&id), sizeof(uint32_t));
        if (!stream.good())
            SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read id of clipped read pair."));
        for (unsigned i = 0; i < clipping.size(); ++i)
        {
            uint8_t tempClip = 0;
            stream.read(reinterpret_cast<char *>(&tempClip), sizeof(unsigned char));
            if (!stream.good())
                SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read number of clipped bases."));
            clipping[i] = tempClip;
        }
        
        clipMap[id] = clipping;
    }
}

// =======================================================================================
// Function readRecords()
// =======================================================================================
// Read all records of one read group in the current window of the input stream. Used by readWindow().
inline void readRecords(seqan::String<ReadPair> & records,
                        zlib_stream::zip_istream & stream,
                        const unsigned & windowBeginPos,
                        const unsigned & numRecords)
{
    for (uint32_t i = 0; i < numRecords; ++i)
    {
        stream.read(reinterpret_cast<char *>(&records[i].pos), sizeof(char));
        if (!stream.good())
            SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read position offset."));
        records[i].pos += windowBeginPos;
        stream.read(reinterpret_cast<char *>(&records[i].distance), sizeof(int32_t));
        if (!stream.good())
            SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read tip distance deviation."));
    }
}

// =======================================================================================
// Function addClippingToRecords()
// =======================================================================================
// Add the clipping info from the map to the records.
inline void addClippingToRecords(seqan::String<ReadPair> & records,
                                 const std::map<unsigned, std::vector<uint8_t>> & clipMap,
                                 const unsigned & numClippedRecords)
{
    if (numClippedRecords == 0)
        return;
    unsigned  processed = 0;
    const std::map<unsigned, std::vector<uint8_t>>::const_iterator mapEnd = clipMap.end();
    for (unsigned  i = 0; i < seqan::length(records) && processed < numClippedRecords; ++i)
    {
        std::map<unsigned, std::vector<uint8_t>>::const_iterator it = clipMap.find(i);
        if (it != mapEnd)
        {
            records[i].clipping = it->second;
            ++processed;
        }
    }
}

// =======================================================================================
// Function readWindow()
// =======================================================================================
// Read the entires of the next window in the opened and unzipped stream.
// Return true on success and false if EOF has been reached.
inline bool readWindow(Window & window,
                       zlib_stream::zip_istream & stream,
                       const unsigned & numReadGroups)
{
    reset(window);
    if (!readCoordinates(window.chrom, window.beginPos, stream))
        return false;

    seqan::resize(window.records, numReadGroups);
    for (unsigned rg = 0; rg < numReadGroups; ++rg)
    {
        std::map<unsigned, std::vector<uint8_t>> clipMap;
        unsigned numRecords = readNumRecords(stream);
        unsigned numClippedRecords = readNumClippedRecords(stream);
        unsigned numNonFRRecords = readNumNonFRRecords(stream);
        seqan::resize(window.records[rg], numRecords);
        readClipping(clipMap, stream, numClippedRecords);
        readOrientations(window.records[rg], stream, numNonFRRecords);
        readRecords(window.records[rg], stream, window.beginPos, numRecords);
        addClippingToRecords(window.records[rg], clipMap, numClippedRecords);
    }
    return true;
}


// -----------------------------------------------------------------------------
// Navigation
// -----------------------------------------------------------------------------


// =======================================================================================
// Function indexBeginPos()
// =======================================================================================
// Return the file position where the index starts.
// FilePos is calculated by adding the bytes of the Magic string, Version, l_region, n_file_offsets
// and n_transloc_file_offsets
inline unsigned indexBeginPos()
{
    return 7 + sizeof(uint16_t) + 3 * sizeof(uint32_t);
}

// =======================================================================================
// Function jumpToRegion()
// =======================================================================================
template<typename TStream>
inline void jumpToRegion(TStream & stream,
                         const seqan::String<seqan::CharString> & contigNames,
                         const seqan::String<int32_t> & contigLengths,
                         unsigned indexRegionSize,
                         seqan::CharString seqName, 
                         uint32_t beginPos)
{
    uint64_t pos = 0;
    stream.clear();

    // Add index positions for other chromosomes before the region.
    for (unsigned i = 0; i < seqan::length(contigNames); ++i)
    {
        SEQAN_ASSERT_LT(i, seqan::length(contigNames));    // Tested before by checkRois().
        if (seqName == contigNames[i])
        {
            break;
        }
        else
        {
            pos += contigLengths[i] / indexRegionSize + 1;
        }
    }
    // Add index positions within the region's chromosome.
    pos += beginPos / indexRegionSize;

    // Go to the position of the region's file offset in index.
    uint64_t offset = 0;
    stream.seekg(indexBeginPos() + pos * sizeof(uint64_t));

    // Read the offset and move the stream.
    stream.read(reinterpret_cast<char *>(&offset), sizeof(uint64_t));
    stream.seekg(offset);
}


//
// Translocation Stuff
//

// =======================================================================================
// Class TranslocationRead
// =======================================================================================
// Holds the info defining a single read of a read pair affected by a translocation
struct TranslocationRead
{
    uint32_t        refID;          // ID of the chromsome the read is mapped to
    uint32_t        pos;            // Position on the reference
    unsigned char   clip_0;           // Number of soft clipped bases at the 5'-end of the read.
    unsigned char   clip_1;           // Number of soft clipped bases at the 3'-end of the read.

    TranslocationRead(): refID(seqan::maxValue<uint32_t>()), pos(0), clip_0(0), clip_1(0) {}

    TranslocationRead(unsigned r, unsigned p, unsigned c_0, unsigned c_1): refID(r), pos(p), clip_0(c_0), clip_1(c_1) {}

    TranslocationRead& operator= (TranslocationRead other)
    {
        std::swap(refID, other.refID);
        std::swap(pos, other.pos);
        std::swap(clip_0, other.clip_0);
        std::swap(clip_1, other.clip_1);
        return *this;
    }
};

// =======================================================================================
// Struct FirstTranslocationRead
// =======================================================================================
// Similar to TranclocationRead, but without member for refID because this info will be stored in the window.
// Should be used in the class TranslocationWindowEntry.
struct FirstTranslocationRead
{
    uint32_t        pos;            // Position on the reference
    unsigned char   clip_0;           // Number of soft clipped bases at the 5'-end of the read.
    unsigned char   clip_1;           // Number of soft clipped bases at the 3'-end of the read.

    FirstTranslocationRead(): pos(0), clip_0(0), clip_1(0) {}
};
// =======================================================================================
// Struct FirstTranslocationRead
// =======================================================================================
// Similar to TranslocationPair, but with a FirstTranslocationRead object for first to avoid storing of 
// redundant refID information in every first reads.
// Should be used in the class TranslocationWindow.
struct TranslocationWindowEntry
{
    FirstTranslocationRead first;
    TranslocationRead second;
    Orientation orientation;            // Orientation of first and second, in this order.

    TranslocationWindowEntry():
        first(), second(), orientation(Orientation::FR)
    {}

    TranslocationWindowEntry(unsigned p, unsigned c, Orientation o)
    {
        first.pos = p;
        first.clip_1 = c;
        orientation = o;
    }
    TranslocationWindowEntry(unsigned p, unsigned c_0, unsigned c_1, Orientation o)
    {
        first.pos = p;
        first.clip_0 = c_0;
        first.clip_1 = c_1;
        orientation = o;
    }
    TranslocationWindowEntry(const FirstTranslocationRead & f, const TranslocationRead & s, Orientation o)
    {
        first = f;
        second = s;
        orientation = o;
    }
};

// =======================================================================================
// Struct TranslocationWindow
// =======================================================================================
// Defines a translocationWindow window on the Genome. Can contain samples from multiple read groups.
struct TranslocationWindow
{
    int32_t chrom;
    int32_t beginPos;
    seqan::String<seqan::String<TranslocationWindowEntry> > records;

    TranslocationWindow():
        chrom(seqan::GenomicRegion::INVALID_ID), beginPos(seqan::GenomicRegion::INVALID_POS)
    {}

    TranslocationWindow(int32_t c, int32_t pos):
        chrom(c), beginPos(pos)
    {}
};

// =======================================================================================
// Function readTranslocationCoordinates()
// =======================================================================================
// Read the chromosome and position of the translocation window from the input stream.
// Return true on success, false if EOF has been reached.
inline bool readTranslocationCoordinates(int32_t & chrom, int32_t & beginPos, zlib_stream::zip_istream & stream)
{
    // Read chromosome and position.
    stream.read(reinterpret_cast<char *>(&chrom), sizeof(uint32_t));
    if (stream.eof())
        return false;
    if (!stream.good())
        SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read contig name in translocation block."));
    stream.read(reinterpret_cast<char *>(&beginPos), sizeof(uint32_t));
    if (!stream.good())
        SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read position in translocation block."));
    return true;
}

// =======================================================================================
// Function unpackOrientations()
// =======================================================================================
// Put the orientations in groups of 4, s.t. they can fill all 8 bits of a uint8_t. Put these uints8_t into a String.
// Overload of unpackOrientations in window_popdel.h for translocation pairs. Does not rely on IDs.
inline void unpackOrientations(seqan::String<TranslocationWindowEntry> & records, const seqan::String<uint8_t> & packed)
{
    uint8_t mask = 192;     // Translates to 11000000 in binary format.
    for (unsigned i = 0; i < length(records); ++i)
        records[i].orientation = static_cast<Orientation>((packed[i / 4] & (mask >> i % 4 * 2)) >> (6-i%4*2));
}

// =======================================================================================
// Function readOrientations()
// =======================================================================================
// Read the packed string of read pair orientations from the stream and unpack them into the records.
// Overload of readOrientations in window.h for translocation windows
inline void readOrientations(seqan::String<TranslocationWindowEntry> & records,
                             zlib_stream::zip_istream & stream,
                             const unsigned numRecords)
{
    if (numRecords == 0)
        return;

    seqan::String<uint8_t> orientations;
    unsigned orientationChars = std::ceil(static_cast<double>(numRecords) / 4);
    seqan::resize(orientations, orientationChars, seqan::Generous()); // TODO: Why Generous()? Try Exact()!
    stream.read(reinterpret_cast<char *>(&orientations[0]), orientationChars);
    if (!stream.good())
        SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read orientations of read pairs in translocation block."));
    unpackOrientations(records, orientations);   // Unpack the orientations and store them in the window.
}
// =======================================================================================
// Function readRecords()
// =======================================================================================
// Read all records of one read group in the current window of the input stream. Used by readWindow().
inline void readRecords(seqan::String<TranslocationWindowEntry> & records,
                        zlib_stream::zip_istream & stream,
                        const unsigned windowBeginPos,
                        const unsigned numRecords)
{
    for (uint32_t i = 0; i < numRecords; ++i)
    {
        stream.read(reinterpret_cast<char *>(&records[i].first.pos), sizeof(char));
        if (!stream.good())
            SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read position offset in translocation block."));
        records[i].first.pos += windowBeginPos;
        stream.read(reinterpret_cast<char *>(&records[i].first.clip_0), sizeof(unsigned char));
        if (!stream.good())
            SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read clipping of first translocation read."));
        stream.read(reinterpret_cast<char *>(&records[i].first.clip_1), sizeof(unsigned char));
        if (!stream.good())
            SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read clipping of first translocation read."));
        stream.read(reinterpret_cast<char *>(&records[i].second.refID), sizeof(uint32_t));
        if (!stream.good())
            SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read refID of second translocation read."));
        stream.read(reinterpret_cast<char *>(&records[i].second.pos), sizeof(uint32_t));
        if (!stream.good())
            SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read position of second translocation read."));
        stream.read(reinterpret_cast<char *>(&records[i].second.clip_0), sizeof(unsigned char));
        if (!stream.good())
            SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read clipping of second translocation read."));
        stream.read(reinterpret_cast<char *>(&records[i].second.clip_1), sizeof(unsigned char));
        if (!stream.good())
            SEQAN_THROW(seqan::ParseError("[PopDel] Unable to read clipping of second translocation read."));
    }
}

// =======================================================================================
// Function readTranslocationWindow()
// =======================================================================================
// Reset the window and get the translocation records from the stream
// Return true on success, false if EOF is encountered
inline bool readTranslocationWindow(TranslocationWindow & window,
                                    zlib_stream::zip_istream & stream,
                                    const unsigned numReadGroups)
{
    reset(window);
    if (!readTranslocationCoordinates(window.chrom, window.beginPos, stream))
        return false;
    seqan::resize(window.records, numReadGroups);
    for (unsigned rg = 0; rg < numReadGroups; ++rg)
    {
        unsigned numRecords = readNumRecords(stream);
        seqan::resize(window.records[rg], numRecords);
        readOrientations(window.records[rg], stream, numRecords);
        readRecords(window.records[rg], stream, window.beginPos, numRecords);
    }
    return true;
}

// =======================================================================================
// Function readTillRoi()
// =======================================================================================
// Read all translocation entries until pos is reached, ignoring everything befor the ROI.
// Return 0 on success, 2 if a new chromsome has been reached before reaching pos and 3 at EOF.
inline unsigned readTillRoi(TranslocationWindow & window,
                            zlib_stream::zip_istream & unzipper,
                            const seqan::CharString & chrom,
                            const unsigned &  pos,
                            const seqan::String<seqan::CharString> & contigNames,
                            int nRG)
{
    do
    {
        if (!readTranslocationWindow(window, unzipper, nRG))
        {
            return 3;               // EOF
        }
        if (contigNames[window.chrom] != chrom)
        {
            return 2;               // New chromosome
        }
    }
    while (static_cast<unsigned>(window.beginPos + 255) < pos);
    return 0;
}

// =======================================================================================
// Function translocationIndexBeginPos()
// =======================================================================================
// Return the file position where the translocation index starts.
// FilePos is calculated by taking the begin pos of the normal index and adding the 
// size of the normal index.
inline unsigned translocationIndexBeginPos(const unsigned numNormalIndexRegions)
{
    return indexBeginPos() + numNormalIndexRegions * sizeof(uint64_t);
}

// =======================================================================================
// Function jumpToTranslocBlock()
// =======================================================================================
// Jump to the start of the translocation block (i.e. the first byte after the translocation guard)
template<typename TStream>
inline void jumpToTranslocBlock(TStream & stream, const unsigned indexOffset)
{
    // Move to the position of the first entry in the translocation index.
    stream.clear();
    
    uint64_t offset = 0;
    stream.seekg(indexOffset);

    // Read the offset and move the stream.
    stream.read(reinterpret_cast<char *>(&offset), sizeof(uint64_t));
    stream.seekg(offset);
}

#endif