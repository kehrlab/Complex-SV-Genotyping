#include "variantMap.hpp"
#include "genomicRegion.hpp"

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

void VariantMap::detectSmallDeletions()
{
    if (this->regions.size() == 0)
        return;
    
    int delSize = -1;
    for (uint32_t i = 0; i + 1 < this->regions.size(); ++i)
        if (isDeletion(delSize, this->regions[i], this->regions[i+1]))
            this->deletionIndices[i] = delSize;
}

bool VariantMap::isDeletion(int & delSize, GenomicRegion & r1, GenomicRegion & r2)
{
    int maxDelSize = 50;
    delSize = -1;
    
    if (r1.getReferenceName() != r2.getReferenceName())
        return false;
    
    
    if (!r1.isReverse() && !r2.isReverse())
    {
        if (r2.getRegionStart() - r1.getRegionEnd() <= maxDelSize && r2.getRegionStart() > r1.getRegionEnd())
        {
            delSize = r2.getRegionStart() - r1.getRegionEnd();
            return true;
        }
    }
    else if (r1.isReverse() && r2.isReverse())
    {
        if (r1.getRegionStart() - r2.getRegionEnd() <= maxDelSize && r1.getRegionStart() > r2.getRegionEnd())
        {
            delSize = r1.getRegionStart() - r2.getRegionEnd();
            return true;
        }
    }
    
    return false;
}

ReadTemplate VariantMap::simulateTemplate(int startIndex, int64_t insertSize, int readLength, ContigInfo & cInfo)
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

    ReadTemplate sT(records, cInfo);
    if (firstPossibleSplit || lastPossibleSplit)
	    sT.markSuspectedSplit();
    return sT;
}

void VariantMap::getGCContentVector(std::vector<float> & gcVector, uint32_t begin, int minInsertSize, int maxInsertSize)
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
        for (uint32_t i = 0; i + 1 < gcVector.size(); ++i)
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
    if (index < this->regionIndices[regionIndex])
        return ("N");
    
    uint32_t relativePosition = index - this->regionIndices[regionIndex];
    if (relativePosition < seqan::length(this->regions[regionIndex].getSequence()))
        return this->regions[regionIndex].getSequence()[relativePosition];
    return seqan::Dna5String("N");
}

std::vector<BamRecord> VariantMap::getReadsFromIndices(int beginIndex, int endIndex, bool isReverse, bool first)
{
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
    
    if (firstIndex == lastIndex || ((std::abs(firstIndex - lastIndex) == 1) && this->deletionIndices.find(firstIndex) != this->deletionIndices.end()))
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

        // consider gapped alignment for small deletions
        bool deletion = false;
        int delSize = -1;
        if (firstIndex != lastIndex) {
            deletion = true;
            delSize = this->deletionIndices[firstIndex];
            this->possibleSplitRead = true;
        }
        BamRecord record(refName, "tempTemplate", start, end, 60, endIndex-beginIndex + 1, clipRight, clipLeft, reverse, true, first, !first, deletion, delSize);
        createdRecords.push_back(record);
    } 
    else if ((std::abs(firstIndex - lastIndex) == 1) && this->deletionIndices.find(firstIndex) == this->deletionIndices.end())
    {
        int startLeft, startRight, endLeft, endRight, clipLeftLeft, clipRightLeft, clipLeftRight, clipRightRight;
        bool reverseLeft, reverseRight;
        bool primaryLeft = false;
        bool primaryRight = false;

        std::string refNameLeft = firstPosition.getReferenceName();
        std::string refNameRight = lastPosition.getReferenceName();

        if (!this->regions[firstIndex].isReverse())
        {
            startLeft = firstPosition.getRegionStart();
            endLeft = this->regions[firstIndex].getRegionEnd();
            clipLeftLeft = 0;
            clipRightLeft = std::abs(lastDistFromBegin) + 1;
            reverseLeft = false;
        } else {
            startLeft = this->regions[firstIndex].getRegionStart();
            endLeft = firstPosition.getRegionStart();
            clipRightLeft = 0;
            clipLeftLeft = std::abs(lastDistFromBegin) + 1;
            reverseLeft = true;
        }
        if (!this->regions[lastIndex].isReverse())
        {
            startRight = this->regions[lastIndex].getRegionStart();
            endRight = lastPosition.getRegionStart();
            clipLeftRight = std::abs(firstDistFromEnd) + 1;
            clipRightRight = 0;
            reverseRight = false;
        } else {
            startRight = lastPosition.getRegionStart();
            endRight = this->regions[lastIndex].getRegionEnd();
            clipLeftRight = 0;
            clipRightRight = std::abs(firstDistFromEnd) + 1;
            reverseRight = true;
        }

        if (std::abs(firstDistFromEnd) >= std::abs(lastDistFromBegin))
            primaryLeft = true;
        else
            primaryRight = true;

        
        if (isReverse)
        {
            reverseLeft = !reverseLeft;
            reverseRight = !reverseRight;
        }
        
	    this->possibleSplitRead = true;
        BamRecord r1(refNameLeft, "tempTemplate", startLeft, endLeft, 60, endIndex - beginIndex + 1, clipRightLeft, clipLeftLeft, reverseLeft, primaryLeft, first, !first);
        BamRecord r2(refNameRight, "tempTemplate", startRight, endRight, 60, endIndex - beginIndex + 1, clipRightRight, clipLeftRight, reverseRight, primaryRight, first, !first);
        createdRecords.push_back(r1);
        createdRecords.push_back(r2);
    } else {
        // find the region that contains most of the read

        // get lengths of read parts
        struct partInfo
        {
            int length;
            int startIndex;
            int lastIndex;
        };
        std::vector<partInfo> parts;

        partInfo temp;
        temp.length = std::abs(firstDistFromEnd);
        temp.startIndex = firstIndex;
        temp.lastIndex = firstIndex;

        for (int i = firstIndex; i < lastIndex; ++i)
        {
            if (this->deletionIndices.find(firstIndex) != this->deletionIndices.end())
            {
                if (i+1 == lastIndex)
                    temp.length += std::abs(lastDistFromBegin);
                else
                    temp.length += this->regionLengths[i+1];
                temp.lastIndex = i+1;
            } else {
                parts.push_back(temp);
                temp.startIndex = i+1;
                temp.lastIndex = i+1;
                if (i+1 == lastIndex)
                    temp.length = std::abs(lastDistFromBegin);
                else
                    temp.length = this->regionLengths[i+1];
            }
        }
        parts.push_back(temp);
        


        // get index of region where the read starts
        uint32_t largestIdx = 0;
        for (uint32_t i = 0; i < parts.size(); ++i)
            if (parts[i].length > parts[largestIdx].length)
                largestIdx = i;


        // extract information for largest part
        std::string refName {""};
        int clipLeft {0}, clipRight {0}, start{-1}, end{-1};
        bool reverse {false};

        refName = this->regions[parts[largestIdx].startIndex].getReferenceName();
        reverse = this->regions[parts[largestIdx].startIndex].isReverse();

        if (parts[largestIdx].startIndex == firstIndex)
        {
            clipLeft = 0;
            if (!reverse)
                start = firstPosition.getRegionStart();
            else
                end = firstPosition.getRegionStart();
        } 
        else 
        {
            clipLeft = 0;
            for (int i = (int) largestIdx - 1; i >= 0; --i)
                clipLeft += parts[i].length;
            if (!reverse)
                start = this->regions[parts[largestIdx].startIndex].getRegionStart();
            else
                end = this->regions[parts[largestIdx].startIndex].getRegionEnd();
        }

        if (parts[largestIdx].lastIndex == lastIndex)
        {
            clipRight = 0;
            if (!reverse)
                end = lastPosition.getRegionStart();
            else
                start = lastPosition.getRegionStart();
        }
        else 
        {
            clipRight = 0;
            for (uint32_t i = largestIdx + 1; i < parts.size(); ++i)
                clipRight += parts[i].length;
            if (!reverse)
                end = this->regions[parts[largestIdx].lastIndex].getRegionEnd();
            else
                start = this->regions[parts[largestIdx].lastIndex].getRegionStart();
        }
        if (isReverse)
        {
            if (reverse)
                reverse = false;
            else
                reverse = true;
        }

        // consider gaps
        bool deletion = false;
        int delSize = -1;
        if (parts[largestIdx].lastIndex - parts[largestIdx].startIndex > 0)
        {
            deletion = true;
            for (int i = parts[largestIdx].startIndex; i < parts[largestIdx].lastIndex; ++i)
                if (this->deletionIndices[i] > 0)
                    delSize = std::max(delSize, this->deletionIndices[i]);
        }

	    this->possibleSplitRead = true;
        BamRecord record(refName, "tempTemplate", start, end, 60, endIndex - beginIndex + 1, clipRight, clipLeft, reverse, true, first, !first, deletion, delSize);
        createdRecords.push_back(record);
    }
    return createdRecords;
}

GenomicRegion VariantMap::getPositionOnReference(int index)
{
    GenomicRegion r;
    int position = 0;
    for (uint32_t i = 0; i + 1 < this->regions.size(); ++i)
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
    int rIdx = ((int) this->regionIndices.size()) - 1;

    if (index >= this->regionIndices[rIdx] && index < (int) this->totalLength)
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
    for (uint32_t i = 0; i + 1 < this->regions.size(); ++i)
        if (position >= this->regionIndices[i] && position < this->regionIndices[i+1])
            return i;
    if (position >= this->regionIndices[this->regionIndices.size() - 1] && position < (int) this->totalLength)
        idx = this->regionIndices.size() - 1;
    return idx;
}
