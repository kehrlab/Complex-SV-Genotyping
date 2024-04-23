#ifndef RECORDHEADER
#define RECORDHEADER

#include <string>
#include <htslib/sam.h>
#include <iostream>

#include "genomicRegion.hpp"

class BamRecord
{
    std::string referenceName;
    std::string templateName;

    int mappingQuality;

    int seqLength;
    int alignmentLength;
    int clipRight;
    int clipLeft;
    int deletionSize;

    int beginPos;
    int endPos;

    bool primary;
    bool reverse;
    bool first;
    bool last;
    bool allProper;
    bool multiple;
    bool unmapped;
    bool qcNoPass;
    bool duplicate;

    bool deletion;

    uint32_t flag;

    void extractRecordInfo(bam1_t *, bam_hdr_t *);

    public:
    BamRecord();
    BamRecord(bam1_t *, bam_hdr_t *);
    BamRecord(std::string, std::string, int, int, bool, bool, bool);
    BamRecord(std::string, std::string, int, int, int, int, int, int, bool, bool, bool, bool);
    BamRecord(std::string, std::string, int, int, int, int, int, int, bool, bool, bool, bool, bool, int);

    std::string getReferenceName();
    std::string getTemplateName();
    int getStartPos();
    int getEndPos();
    int getFivePrimePos();
    int getThreePrimePos();
    int getSeqLength();
    int getClipRight();
    int getClipLeft();
    int getAlignmentLength();
    bool isPrimary();
    bool isReverse();
    bool isFirst();
    bool isLast();
    bool isClipped();
    bool isClipped(int);
    int getMapQ();
    bool hasFlagAllProper();
    bool containsDeletion();
    int getDeletionSize();

    GenomicRegion getAlignmentRegion();
    bool passesInsertFilter();
    bool passesStandardFilter();
    void print();
};

#endif
