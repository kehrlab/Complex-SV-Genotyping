#include "record.hpp"

BamRecord::BamRecord()
{}

BamRecord::BamRecord(std::string refName, std::string templateName, int begin, int end, bool reverse, bool first, bool last)
{
    this->referenceName = refName;
    this->templateName = templateName;
    this->beginPos = begin;
    this->endPos = end;
    this->reverse = reverse;

    this->mappingQuality = 60;
    this->seqLength = end-begin;
    this->clipLeft = 0;
    this->clipRight = 0;
    this->alignmentLength = this->seqLength;
    this->deletionSize = -1;

    this->primary = true;
    this->first = first;
    this->last = last;
    this->allProper = true;
    this->deletion = false;
}

BamRecord::BamRecord(std::string refName, std::string templateName, int begin, int end, int mappingQuality, int seqLength, int clipRight, int clipLeft, bool reverse, bool primary, bool first, bool last)
{
    this->referenceName = refName;
    this->templateName = templateName;
    this->beginPos = begin;
    this->endPos = end;
    this->alignmentLength = end-begin;
    this->reverse = reverse;
    this->mappingQuality = mappingQuality;
    this->seqLength = seqLength;
    this->clipLeft = clipLeft;
    this->clipRight = clipRight;
    this->primary = primary;
    this->first = first;
    this->last = last;
    this->allProper = true;
    this->deletion = false;
    this->deletionSize = -1;
}

BamRecord::BamRecord(std::string refName, std::string templateName, int begin, int end, int mappingQuality, int seqLength, int clipRight, int clipLeft, bool reverse, bool primary, bool first, bool last, bool deletion, int deletionSize)
{
    this->referenceName = refName;
    this->templateName = templateName;
    this->beginPos = begin;
    this->endPos = end;
    this->alignmentLength = end-begin;
    this->reverse = reverse;
    this->mappingQuality = mappingQuality;
    this->seqLength = seqLength;
    this->clipLeft = clipLeft;
    this->clipRight = clipRight;
    this->primary = primary;
    this->first = first;
    this->last = last;
    this->allProper = true;
    this->deletion = deletion;
    this->deletionSize = deletionSize;
}

BamRecord::BamRecord(bam1_t * record, bam_hdr_t * hdr)
{
    extractRecordInfo(record, hdr);
}

void BamRecord::extractRecordInfo(bam1_t * record, bam_hdr_t * hdr)
{
    this->referenceName = hdr->target_name[record->core.tid];
    this->templateName = bam_get_qname(record);
    this->alignmentLength = bam_cigar2rlen(record->core.n_cigar, bam_get_cigar(record));
    this->beginPos = record->core.pos;
    this->endPos = this->beginPos + this->alignmentLength - 1;
    this->seqLength = record->core.l_qseq;
    this->mappingQuality = record->core.qual;

    // interpret flags
    this->flag = record->core.flag;
    this->first = (bool) (this->flag & BAM_FREAD1);
    this->last = (bool) (this->flag & BAM_FREAD2);
    this->primary = ! ((bool) (this->flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)));
    this->reverse = bam_is_rev(record);
    this->allProper = (bool) (this->flag & BAM_FPROPER_PAIR);
    this->multiple = (bool) (this->flag & BAM_FPAIRED);
    this->duplicate = (bool) (this->flag & BAM_FDUP);
    this->qcNoPass = (bool) (this->flag & BAM_FQCFAIL);
    this->unmapped = (bool) (this->flag & BAM_FMUNMAP);

    // extract CIGAR operations and lenghts
    uint32_t * cigarString = bam_get_cigar(record);
    uint32_t cigar_len = record->core.n_cigar;

    this->clipLeft = 0;
    this->clipRight = 0;
    this->deletion = false;
    this->deletionSize = 0;
    
    for (uint32_t i = 0; i < cigar_len; ++i)
    {
        int oplen = bam_cigar_oplen(cigarString[i]);
        if (((cigarString[i] & 15) == BAM_CSOFT_CLIP) || ((cigarString[i] & 15) == BAM_CHARD_CLIP))
        {
            if (i == 0)
                this->clipLeft = oplen;
            if (i + 1 == cigar_len)
                this->clipRight = oplen;
            continue;
        }
        if ((cigarString[i] & 15) == BAM_CDEL)
        {
            this->deletion = true;
            this->deletionSize = oplen;
            continue;
        }
    }
}

bool BamRecord::hasFlagAllProper()
{
	return this->allProper;
}

std::string BamRecord::getReferenceName()
{
    return this->referenceName;
}

std::string BamRecord::getTemplateName()
{
    return this->templateName;
}

int BamRecord::getStartPos()
{
    return this->beginPos;
}

int BamRecord::getEndPos()
{
    return this->endPos;
}

int BamRecord::getFivePrimePos()
{
    if (!isReverse())
        return getStartPos();
    else
        return getEndPos();
}

int BamRecord::getThreePrimePos()
{
    if (isReverse())
        return getStartPos();
    else
        return getEndPos();
}

int BamRecord::getSeqLength()
{
    return this->seqLength;
}

int BamRecord::getClipLeft()
{
    return this->clipLeft;
}

int BamRecord::getClipRight()
{
    return this->clipRight;
}

int BamRecord::getAlignmentLength()
{
    return this->alignmentLength;
}

bool BamRecord::isPrimary()
{
    return this->primary;
}

bool BamRecord::isReverse()
{
    return this->reverse;
}

bool BamRecord::isFirst()
{
    return this->first && !this->last;
}

bool BamRecord::isLast()
{
    return this->last && !this->first;
}

GenomicRegion BamRecord::getAlignmentRegion()
{
    return GenomicRegion(getReferenceName(), std::min(getFivePrimePos(), getThreePrimePos()), std::max(getFivePrimePos(), getThreePrimePos()));
}

int BamRecord::getMapQ()
{
    return this->mappingQuality;
}

bool BamRecord::isClipped()
{
    return !(this->clipLeft == 0 && this->clipRight == 0);
}

bool BamRecord::containsDeletion()
{
    return this->deletion;
}

int BamRecord::getDeletionSize()
{
    return this->deletionSize;
}

bool BamRecord::passesInsertFilter()
{
    if (this->unmapped || this->qcNoPass || this->duplicate || !this->multiple || !this->primary)
        return false;
    if (isClipped() || containsDeletion())
        return false;
    if (this->mappingQuality < 40)
        return false;
    return true;
}

bool BamRecord::passesStandardFilter()
{
    return !(this->unmapped || this->duplicate || !this->multiple || !this->primary);
}
