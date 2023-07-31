#include "variantMap.hpp"
#include "genomicRegion.hpp"
#include "readTemplate.hpp"
#include "record.hpp"

VariantMap::VariantMap()
{
    this->totalLength = 0;
}

void VariantMap::calculateMapLength()
{
    for (GenomicRegion region : this->regions)
    {
        this->regionIndices.push_back(totalLength);
        int l = std::abs(region.getRegionEnd() - region.getRegionStart()) + 1;
        this->regionLengths.push_back(l);
        this->totalLength += l;
    }
}

ReadTemplate VariantMap::simulateTemplate(int startIndex, int insertSize, int readLength)
{
    std::vector<BamRecord> records;

    std::string templateName = "tempTemplate";

    int beginFirstIndex = startIndex;
    int endFirstIndex;
    if (insertSize >= readLength)
        endFirstIndex = startIndex + readLength - 1;
    else
        endFirstIndex = startIndex + insertSize - 1;
    this->possibleSplitRead = false;
    std::vector<BamRecord> firstRecords = getReadsFromIndices(beginFirstIndex, endFirstIndex, false, true);
    bool firstPossibleSplit = this->possibleSplitRead;

    int beginLastIndex = startIndex + insertSize - readLength;
    if (insertSize < readLength)
        beginLastIndex = startIndex;
    int endLastIndex = startIndex + insertSize - 1;

    this->possibleSplitRead = false;
    std::vector<BamRecord> lastRecords = getReadsFromIndices(beginLastIndex, endLastIndex, true, false);
    bool lastPossibleSplit = this->possibleSplitRead;

    for (BamRecord r : firstRecords)
        records.push_back(r);
    for (BamRecord r : lastRecords)
        records.push_back(r);

    ReadTemplate sT(records);
    if (firstPossibleSplit || lastPossibleSplit)
	    sT.markSuspectedSplit();
    return sT;
}

void VariantMap::getGCContentVector(std::vector<float> & gcVector, unsigned begin, int minInsertSize, int maxInsertSize)
{
    seqan::Dna5String bg = "G";
    seqan::Dna5String bc = "C";
    seqan::Dna5String base;
    int total = 0;
    int gcCount = 0;
    float temp = 0;
    
    if (gcVector.size() == 0)
    {
        for (int i = 0; i <= maxInsertSize; ++i)
        {
            ++total;
            base = getBaseAtIndex(begin + i);
            if (base == bg || base == bc)
                ++gcCount;
            if (i >= minInsertSize)
                gcVector.push_back((float) gcCount / (float) total);
        }
    } else {
        base = getBaseAtIndex(begin - 1);
        total = minInsertSize + 1;
        for (int i = 0; i < gcVector.size() - 1; ++i)
        {
            ++total;
            temp = gcVector[i + 1] * total;
            if (base == bg || base == bc)
                temp -= 1;
            gcVector[i] = temp / (total - 1);
        }
        total = maxInsertSize + 1;
        temp = gcVector[gcVector.size() - 1] * total;
        if (base == bg || base == bc)
            temp -= 1;
        base = getBaseAtIndex(begin + maxInsertSize);
        if (base == bg || base == bc)
            temp += 1;
        gcVector[gcVector.size() - 1] = temp / total;
    }
}

seqan::Dna5String VariantMap::getBaseAtIndex(int index)
{
    int regionIndex = getRegionIndex(index);
    if (regionIndex == -1)
        return seqan::Dna5String("N");
    int relativePosition = index - this->regionIndices[regionIndex];
    if (relativePosition < seqan::length(this->regions[regionIndex].getSequence()))
        return this->regions[regionIndex].getSequence()[relativePosition];
    return seqan::Dna5String("N");
}

std::vector<BamRecord> VariantMap::getReadsFromIndices(int beginIndex, int endIndex, bool isReverse, bool first)
{
    // is there a distinction between clipping at 5' and 3' ends?
    // not sure...
    // there might be a bias towards keeping 5' ends due to higher quality base calls...
    // I think I am still observing edge effects for small complex variants

    std::vector<BamRecord> createdRecords;
    int firstIndex = getRegionIndex(beginIndex);
    int lastIndex = getRegionIndex(endIndex);

    GenomicRegion firstPosition, lastPosition;
    firstPosition = getPositionOnReference(beginIndex);
    lastPosition = getPositionOnReference(endIndex);

    int firstDistFromEnd, lastDistFromBegin;
    if (this->regions[firstIndex].isReverse())
        firstDistFromEnd = firstPosition.getRegionStart() - this->regions[firstIndex].getRegionStart();
    else
        firstDistFromEnd = this->regions[firstIndex].getRegionEnd() - firstPosition.getRegionStart();
    if (this->regions[lastIndex].isReverse())
        lastDistFromBegin = this->regions[lastIndex].getRegionEnd() - lastPosition.getRegionStart();
    else
        lastDistFromBegin = lastPosition.getRegionStart() - this->regions[lastIndex].getRegionStart();
    
    if (firstIndex == lastIndex)
    {
        std::string refName = firstPosition.getReferenceName();
        int start, end, clipLeft, clipRight;
        bool reverse;
        if (this->regions[firstIndex].isReverse())
        {
            start = lastPosition.getRegionStart();
            end = firstPosition.getRegionStart();
            reverse = true;
        } else {
            start = firstPosition.getRegionStart();
            end = lastPosition.getRegionStart();
            reverse = false;
        }
        clipLeft = 0; // no clipping
        clipRight = 0;
        if (isReverse)
        {
            if (reverse)
                reverse = false;
            else
                reverse = true;
        }
        BamRecord record(refName, "tempTemplate", start, end, 60, endIndex-beginIndex, clipRight, clipLeft, reverse, true, first, !first);
        createdRecords.push_back(record);
    } 
    else if (std::abs(firstIndex - lastIndex) == 1)
    {
        std::string refNameLeft, refNameRight;
        int startLeft, startRight, endLeft, endRight, clipLeftLeft, clipRightLeft, clipLeftRight, clipRightRight;
        bool reverseLeft, reverseRight;
        bool primaryLeft = false;
        bool primaryRight = false;

        refNameLeft = firstPosition.getReferenceName();
        refNameRight = lastPosition.getReferenceName();

        if (!this->regions[firstIndex].isReverse())
        {
            startLeft = firstPosition.getRegionStart();
            endLeft = this->regions[firstIndex].getRegionEnd();
            clipLeftLeft = 0;
            clipRightLeft = std::abs(lastDistFromBegin);
            reverseLeft = false;
        } else {
            startLeft = this->regions[firstIndex].getRegionStart();
            endLeft = firstPosition.getRegionStart();
            clipRightLeft = 0;
            clipLeftLeft = std::abs(lastDistFromBegin);
            reverseLeft = true;
        }
        if (!this->regions[lastIndex].isReverse())
        {
            startRight = this->regions[lastIndex].getRegionStart();
            endRight = lastPosition.getRegionStart();
            clipLeftRight = std::abs(firstDistFromEnd);
            clipRightRight = 0;
            reverseRight = false;
        } else {
            startRight = lastPosition.getRegionStart();
            endRight = this->regions[lastIndex].getRegionEnd();
            clipLeftRight = 0;
            clipRightRight = std::abs(firstDistFromEnd);
            reverseRight = true;
        }

        if (std::abs(firstDistFromEnd) >= std::abs(lastDistFromBegin))
            primaryLeft = true;
        else
            primaryRight = true;

        if (isReverse)
        {
            if (reverseLeft)
                reverseLeft = false;
            else
                reverseLeft = true;
            if (reverseRight)
                reverseRight = false;
            else
                reverseRight = true;
        }
        
	this->possibleSplitRead = true;
        BamRecord r1(refNameLeft, "tempTemplate", startLeft, endLeft, 60, endIndex - beginIndex, clipRightLeft, clipLeftLeft, reverseLeft, primaryLeft, first, !first);
        BamRecord r2(refNameRight, "tempTemplate", startRight, endRight, 60, endIndex - beginIndex, clipRightRight, clipLeftRight, reverseRight, primaryRight, first, !first);
        createdRecords.push_back(r1);
        createdRecords.push_back(r2);
    } else {
        std::vector<int> partLengths;
        // find the region that contains most of the read

        // first region
        partLengths.push_back(std::abs(firstDistFromEnd));
        // in between
        for (int i = firstIndex + 1; i < lastIndex; ++i)
        {
            partLengths.push_back(this->regionLengths[i]);
        }
        // last region
        partLengths.push_back(std::abs(lastDistFromBegin));

        int largestIdx = 0;
        for (int i = 0; i < partLengths.size(); ++i)
            if (partLengths[i] > partLengths[largestIdx])
                largestIdx = i;
        int idx = firstIndex + largestIdx; // region index

        // extract information for largest part
        std::string refName;
        int clipLeft, clipRight, start, end;
        bool reverse;

        if (largestIdx == 0)
        {
            // right clipping, left ist the new read
            refName = firstPosition.getReferenceName();
            clipLeft = 0;
            clipRight = 0;
            if (! this->regions[idx].isReverse())
            {
                start = firstPosition.getRegionStart();
                end = this->regions[idx].getRegionEnd();
                for (int i = 1; i < partLengths.size(); ++i)
                    clipRight += partLengths[i];
                reverse = false;
            } else {
                start = this->regions[idx].getRegionStart();
                end = firstPosition.getRegionStart();
                reverse = true;
                for (int i = 1; i < partLengths.size(); ++i)
                    clipLeft += partLengths[i];
            }  
        } else if (largestIdx == partLengths.size() - 1)
        {
            refName = lastPosition.getReferenceName();
            clipRight = 0;
            clipLeft = 0;
            if (! this->regions[idx].isReverse())
            {
                start = this->regions[idx].getRegionStart();
                end = lastPosition.getRegionStart();
                for (int i = 0; i < partLengths.size() - 1; ++i)
                    clipLeft += partLengths[i];
                reverse = false;
            } else {
                start = lastPosition.getRegionStart();
                end = this->regions[idx].getRegionEnd();
                reverse = true;
                for (int i = 0; i < partLengths.size() - 1; ++i)
                    clipRight += partLengths[i];
            }
        } else {
            refName = this->regions[idx].getReferenceName();
            start = this->regions[idx].getRegionStart();
            end = this->regions[idx].getRegionEnd();
            reverse = this->regions[idx].isReverse();
            clipLeft = 0;
            clipRight = 0;
            for (int i = 0; i < partLengths.size(); ++i)
            {
                if (i < largestIdx)
                    clipLeft += partLengths[i];
                if (i > largestIdx)
                    clipRight += partLengths[i];
            }
            if (reverse)
            {
                int temp = clipRight;
                clipRight = clipLeft;
                clipLeft = temp;
            }
        }
        if (isReverse)
        {
            if (reverse)
                reverse = false;
            else
                reverse = true;
        }
	this->possibleSplitRead = true;
        BamRecord record(refName, "tempTemplate", start, end, 60, lastIndex - beginIndex, clipRight, clipLeft, reverse, true, first, !first);
        createdRecords.push_back(record);
    }
    return createdRecords;
}

GenomicRegion VariantMap::getPositionOnReference(int index)
{
    GenomicRegion r;
    int position = 0;
    for (unsigned i = 0; i < this->regions.size() - 1; ++i)
    {
        if (index < this->regionIndices[i] || index >= this->regionIndices[i+1])
            continue;
        
        int remainder = index - this->regionIndices[i];
        
        if (! this->regions[i].isReverse())
            position = this->regions[i].getRegionStart() + remainder;
        else
            position = this->regions[i].getRegionEnd() - remainder;
        
        r.setRegionStart(position);
        r.setRegionEnd(position);
        r.setReferenceName(this->regions[i].getReferenceName());
        return r;
    }
    int rIdx = this->regionIndices.size() - 1;
    if (index >= this->regionIndices[rIdx] && index < this->totalLength)
    {
        int remainder = index - this->regionIndices[rIdx];
        if (! this->regions[rIdx].isReverse())
            position = this->regions[rIdx].getRegionStart() + remainder;
        else
            position = this->regions[rIdx].getRegionEnd() - remainder;
        r.setRegionStart(position);
        r.setRegionEnd(position);
        r.setReferenceName(this->regions[rIdx].getReferenceName());
    }
    return r;
}

int VariantMap::getRegionIndex(int position)
{   
    int idx = -1;
    for (unsigned i = 0; i < this->regions.size() - 1; ++i)
        if (position >= this->regionIndices[i] && position < this->regionIndices[i+1])
            return i;
    if (position >= this->regionIndices[this->regionIndices.size() - 1] && position < this->totalLength)
        idx = this->regionIndices.size() - 1;
    return idx;
}
