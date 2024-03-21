#include "readTemplate.hpp"
#include <algorithm>
#include <stdexcept>

ReadTemplate::ReadTemplate()
{}

ReadTemplate::ReadTemplate(std::vector<BamRecord> records, ContigInfo & cInfo)
{
    this->records = records;
    this->split = false;
    this->spanning = false;
    this->insertSize = 0;
    this->interChromosome = false;
    this->templateWeight = 0;
    this->templateGCContent = 0;
    this->templateName = records[0].getTemplateName();
    this->junctionString = "";
    this->bpSpanString = "";
    this->splitFirst = false;
    this->splitLast = false;
    this->suspectedSplit = false;

    determineInformativeRecords();
    determineFivePrimeEnds(cInfo);

    determineOrientation();
    calculateInsertSize();
}

void ReadTemplate::determineInformativeRecords()
{
    this->primaryFirst = -1;
    this->primaryLast = -1;
    for (uint32_t i = 0; i < this->records.size(); ++i)
    {
        if (records[i].isFirst())
            this->firstIndices.push_back((int) i);
        else if (records[i].isLast())
            this->lastIndices.push_back((int) i);
    }

    for (auto index : this->firstIndices)
        if (this->records[index].isPrimary())
            this->primaryFirst = index;
    for (auto index : this->lastIndices)
        if (this->records[index].isPrimary())
            this->primaryLast = index;
    if (isProperPair())
    {
        float q1 = std::pow(10, - this->records[this->primaryFirst].getMapQ()/10);
        float q2 = std::pow(10, - this->records[this->primaryLast].getMapQ()/10);
        this->templateWeight = (1-q1)*(1-q2);
    }
}

void ReadTemplate::determineFivePrimeEnds(ContigInfo & cInfo)
{
    if (this->primaryFirst == -1 || this->primaryLast == -1)
    {
        this->fivePrimeFirst = -1;
        this->fivePrimeLast = -1;
        return;
    }

    this->fivePrimeFirst = this->records[this->primaryFirst].getFivePrimePos();
    if (this->records[this->primaryFirst].isReverse())
        this->fivePrimeFirst += this->records[this->primaryFirst].getClipRight();
    else
        this->fivePrimeFirst -= this->records[this->primaryFirst].getClipLeft();

    this->fivePrimeLast = this->records[this->primaryLast].getFivePrimePos();
    if (this->records[this->primaryLast].isReverse())
        this->fivePrimeLast += this->records[this->primaryLast].getClipRight();
    else
        this->fivePrimeLast -= this->records[this->primaryLast].getClipLeft();

    // switch 5' positions to global coordinates in order to define insert size for inter-chromosomal read pairs
    if (isProperPair()) {
        if (cInfo.globalPositions.find(this->records[this->primaryFirst].getReferenceName()) == cInfo.globalPositions.end() || cInfo.globalPositions.find(this->records[this->primaryLast].getReferenceName()) == cInfo.globalPositions.end())
        {
            std::cerr << "One of the chromosomes was not contained in cInfo" << std::endl;
            std::cerr << this->records[this->primaryFirst].getReferenceName() << std::endl;
            std::cerr << this->records[this->primaryLast].getReferenceName() << std::endl << std::endl;
            std::cerr << "Present contigs: " << std::endl;
            for (auto c : cInfo.globalPositions)
                std::cerr << c.first << ":\t" << c.second << std::endl;
            std::cerr << "-----------------------------------" << std::endl;
            
            if (this->records[this->primaryFirst].getReferenceName() == this->records[this->primaryLast].getReferenceName())
                throw std::runtime_error("Could not calculate insert size for inter-chromosomal read pair. Exit.");
        } else {
            this->fivePrimeFirst = cInfo.globalPositions[this->records[this->primaryFirst].getReferenceName()] + this->fivePrimeFirst;
            this->fivePrimeLast = cInfo.globalPositions[this->records[this->primaryLast].getReferenceName()] + this->fivePrimeLast;
        }
    }
}

void ReadTemplate::determineOrientation()
{
    if (this->primaryFirst == -1 || this->primaryLast == -1)
    {
        this->orientation = "";
        return;
    }

    this->orientation = "FR";

    if (this->records[this->primaryFirst].isReverse() && this->records[this->primaryLast].isReverse())
        this->orientation = "RR";
    else if (!this->records[this->primaryFirst].isReverse() && !this->records[this->primaryLast].isReverse())
        this->orientation = "FF";
}

void ReadTemplate::calculateInsertSize()
{
    this->insertSize = 0;

    if (isProperPair())
    {   
        if (this->orientation == "RR" || this->orientation == "FF")
        {
            this->insertSize = std::abs(this->fivePrimeFirst - this->fivePrimeLast) + 1;
            return;
        }
        
        if (this->orientation == "FR" || this->orientation == "RF")
        {   
            if (this->records[this->primaryFirst].isReverse() && !this->records[this->primaryLast].isReverse())
                this->insertSize = this->fivePrimeFirst - this->fivePrimeLast + 1;
            else if (!this->records[this->primaryFirst].isReverse() && this->records[this->primaryLast].isReverse())
                this->insertSize = this->fivePrimeLast - this->fivePrimeFirst + 1;
            return;
        }
    }
}

void ReadTemplate::determineLocationStrings()
{
    
    if (this->splittingJunctions.size() > 0)
    {
        std::vector<int> junctions;
        for (auto id : this->splittingJunctions)
            junctions.push_back(id);
        std::sort(junctions.begin(), junctions.end());
        this->junctionString = "";
        for (int i : junctions)
            this->junctionString = this->junctionString + std::to_string(i);
    }
    
    if (this->spanningBreakpoints.size() > 0)
    {
        std::vector<int> breakpoints;
        for (auto id : this->spanningBreakpoints)
            breakpoints.push_back(id);
        std::sort(breakpoints.begin(), breakpoints.end());
        this->bpSpanString = "";
        for (int i : breakpoints)
            this->bpSpanString = this->bpSpanString + std::to_string(i);
    }

    if (this->bridgedBreakpoints.size() > 0)
    {
        std::vector<int> breakpoints;
        for (auto id : this->bridgedBreakpoints)
            breakpoints.push_back(id);
        std::sort(breakpoints.begin(), breakpoints.end());
        this->bpBridgeString = "";
        for (int i : breakpoints)
            this->bpBridgeString = this->bpBridgeString + std::to_string(i);
    }
}

void ReadTemplate::findSplitReads(std::vector<Junction> & junctions, std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> & chromosomeStructures, int sMax)
{
    // both reads must properly align
    if (!isProperPair())
        return;
    
    // there are two modes of detection / two different kinds of splits
    findSplitsBasedOnClipping(chromosomeStructures, sMax);
    findSplitsBasedOnGaps(junctions);
    return;
}


void ReadTemplate::findSplitsBasedOnClipping(std::unordered_map<std::string, std::unordered_map<std::string, JunctionRegion>> & chromosomeStructures, int sMax)
{
    BamRecord rFirst = this->records[this->primaryFirst];
    BamRecord rLast = this->records[this->primaryLast];

    // only works with at least one clipped record 
    if (!rFirst.isClipped() && !rLast.isClipped())
	    return;

    SplitAlignmentInfo splitInfo;

    // go over all chromosomes on all alleles
    for (auto & alleleStruct : chromosomeStructures)
        if (alleleStruct.first != "REF")
            for (auto chr : alleleStruct.second)
                findClippedSplitsOnChromosome(chr.first, chr.second, splitInfo, sMax);
    
    // evaluate the results
    if (splitInfo.junctionIndices.size() > 1) {
        this->junctionString = "ambiguous";
        this->split = true;
    } else if (splitInfo.junctionIndices.size() == 1) {
        for (auto & idx : splitInfo.junctionIndices[0])
            this->splittingJunctions.insert(idx);
        this->split = true;
    }

    return;    
}

void ReadTemplate::findClippedSplitsOnChromosome(std::string cName, JunctionRegion & jRegion, SplitAlignmentInfo & splitInfo, int sMax)
{
    BamRecord & rFirst = this->records[this->primaryFirst];
    BamRecord & rLast = this->records[this->primaryLast];

    // find possible locations and associated junctions for the two records
    std::vector<GenomicRegion> firstRegions;
    std::vector<std::vector<int>> firstJunctionIndices; 
    getIndexRegions(firstRegions, firstJunctionIndices, rFirst, jRegion);

    std::vector<GenomicRegion> lastRegions;
    std::vector<std::vector<int>> lastJunctionIndices; 
    getIndexRegions(lastRegions, lastJunctionIndices, rLast, jRegion);

    std::vector<std::unordered_set<int>> indexGroups;
    std::vector<int> insertSizes;

    // determine possible combinations and their orientation and insert size on the variant allele(s)
    for (uint32_t i = 0; i < firstRegions.size(); ++i)
    {
        for (uint32_t j = 0; j < lastRegions.size(); ++j)
        {
            std::unordered_set<int> indices;
            int64_t s = 0;

            // skip read pairs with invalid original orientations
            if (firstRegions[i].isReverse() == lastRegions[j].isReverse())
                continue;
            
            if (firstRegions[i].isReverse() && !lastRegions[j].isReverse())
            {
                // skip read pair if original orientation is RF
                if (firstRegions[i].getRegionStart() < lastRegions[j].getRegionStart())
                    continue;
                
                // calculate insert size on variant allele
                s = firstRegions[i].getRegionEnd() - lastRegions[j].getRegionStart() + 1;
                // if it is not within the range of the insert size distribution, skip
                if (s < 0 || s > sMax)
                    continue;

                for (auto & idx : firstJunctionIndices[i])
                    indices.insert(idx);
                for (auto & idx : lastJunctionIndices[j])
                    indices.insert(idx);
            }

            if (!firstRegions[i].isReverse() && lastRegions[j].isReverse())
            {
                // skip read pair if original orientation is RF
                if (firstRegions[i].getRegionStart() > lastRegions[j].getRegionStart())
                    continue;
                
                // calculate insert size
                s = lastRegions[j].getRegionEnd() - firstRegions[i].getRegionStart() + 1;
                // if it is not within the range of the insert size distribution, skip
                if (s < 0 || s > sMax)
                    continue;

                for (auto & idx : firstJunctionIndices[i])
                    indices.insert(idx);
                for (auto & idx : lastJunctionIndices[j])
                    indices.insert(idx);
            }
            
            if (indices.size() > 0)
            {	
                indexGroups.push_back(indices);
                insertSizes.push_back(s);
            }
	    }	
    }
    // store all combinations and deal with them when adding to profile
    splitInfo.junctionIndices = indexGroups;
    splitInfo.insertSize = insertSizes;
}


void ReadTemplate::getIndexRegions(std::vector<GenomicRegion> & regions, std::vector<std::vector<int>> & indices, BamRecord & record, JunctionRegion & jRegion)
{
    int overlap = 20;
    int smallMargin = 10;
    GenomicRegion recordRegion {record.getAlignmentRegion()};

    // check for each region whether read overlaps
    for (uint32_t i = 0; i < jRegion.regions.size(); ++i)
    {
        GenomicRegion & r{jRegion.regions[i]};

        if (!r.overlaps(recordRegion, smallMargin))
            continue;

        int beginIdx {-1}, endIdx {-1};
        std::vector<int> tempIndices;

        // check for clipping at the end of first region
        if (i == 0)
        {
            beginIdx = jRegion.junctionIndices[0];
            endIdx = jRegion.junctionIndices[0];

            int clipRight{0};
            int dRight{0};

            if (!r.isReverse()) {
                dRight = r.getRegionEnd() - recordRegion.getRegionEnd();
                beginIdx -= (r.getRegionEnd() - recordRegion.getRegionStart());
                endIdx -= dRight;
                clipRight = record.getClipRight();
            } else {
                dRight = recordRegion.getRegionStart() - r.getRegionStart();
                beginIdx -= (recordRegion.getRegionEnd() - r.getRegionStart());
                endIdx -= dRight;
                clipRight = record.getClipLeft();
            }

            if (clipRight > 0 && dRight >= 0 && clipRight - dRight >= overlap && dRight < smallMargin)
            {
                endIdx += clipRight;
                tempIndices.push_back(jRegion.junctions[i].getID());
            }
        }
        
        // check for clipping at each end of a middle region
        if (i > 0 && i + 1 < jRegion.regions.size())
        {
            beginIdx = jRegion.junctionIndices[i-1];
            endIdx = jRegion.junctionIndices[i-1];

            int clipLeft{0}, clipRight{0};
            int dLeft{0}, dRight{0};

            if (!r.isReverse()) {
                dLeft = recordRegion.getRegionStart() - r.getRegionStart();
                dRight = r.getRegionEnd() - recordRegion.getRegionEnd();
                beginIdx += dLeft;
                endIdx += (recordRegion.getRegionEnd() - r.getRegionStart());
                clipLeft = record.getClipLeft();
                clipRight = record.getClipRight();
            } else {
                dLeft = r.getRegionEnd() - recordRegion.getRegionEnd();
                dRight = recordRegion.getRegionStart() - r.getRegionStart();
                beginIdx += dLeft;
                endIdx += (r.getRegionEnd() - recordRegion.getRegionStart());
                clipRight = record.getClipLeft();
                clipLeft = record.getClipRight();
            }

            
            if (clipLeft > 0 && dLeft >= 0 && clipLeft - dLeft >= overlap && dLeft < smallMargin)
            {
                beginIdx -= clipLeft;
                tempIndices.push_back(jRegion.junctions[i-1].getID());
            }

            if (clipRight > 0 && dRight >= 0 && clipRight - dRight >= overlap && dRight < smallMargin)
            {
                endIdx += clipRight;
                tempIndices.push_back(jRegion.junctions[i].getID());
            }
        }

        // check for clipping on the left side of last region
        uint32_t nJ = jRegion.regions.size() - 1;
        if (i == nJ)
        {
            beginIdx = jRegion.junctionIndices[nJ-1];
            endIdx = jRegion.junctionIndices[nJ-1];

            int clipLeft{0};
            int dLeft{0};

            if (!r.isReverse()) {
                dLeft = recordRegion.getRegionStart() - r.getRegionStart();
                beginIdx += dLeft;
                endIdx += recordRegion.getRegionEnd() - r.getRegionStart();
                clipLeft = record.getClipLeft();
            } else {
                dLeft = r.getRegionEnd() - recordRegion.getRegionEnd();
                beginIdx += dLeft;
                endIdx += r.getRegionEnd() - recordRegion.getRegionStart();
                clipLeft = record.getClipRight();
            }

            if (clipLeft > 0 && dLeft >= 0 && clipLeft - dLeft >= overlap && dLeft < smallMargin)
            {
                beginIdx -= clipLeft;
                tempIndices.push_back(jRegion.junctions[nJ-1].getID());
            }
        }

        // infer original read orientation
        bool originalReverse = (r.isReverse() != record.isReverse());
        
        regions.push_back(GenomicRegion{jRegion.chromosome, beginIdx, endIdx, originalReverse});
        indices.push_back(tempIndices);
    }
}


bool ReadTemplate::alignsWithinExpectedDistance(std::vector<Junction> & junctions, uint32_t index, BamRecord & clippedRecord, BamRecord & secondRecord)
{
    int expectedDistance = 1500;
    int remainingDistance = expectedDistance;

    // left 
    for (int i = (int) index; i > 0; --i)
    {
        GenomicRegion junctionRegion(
            junctions[i].getRefNameLeft(),
            std::min(junctions[i].getPositionLeft(), junctions[i-1].getPositionRight()),
            std::max(junctions[i].getPositionLeft(), junctions[i-1].getPositionRight())
        );
        if (alignsWithinRegion(remainingDistance, junctionRegion, secondRecord, junctions[i].getPositionLeft()))
            return true;
        if (remainingDistance <= 0)
            break;
    }
    if (remainingDistance > 0)
    {
        GenomicRegion junctionRegion(
            junctions[0].getRefNameLeft(),
            std::min(junctions[0].getPositionLeft(), junctions[0].getPositionLeft() + junctions[0].getDirectionLeft() * 1000),
            std::max(junctions[0].getPositionLeft(), junctions[0].getPositionLeft() + junctions[0].getDirectionLeft() * 1000)
        );
        if (alignsWithinRegion(remainingDistance, junctionRegion, secondRecord, junctions[0].getPositionLeft()))
            return true;
    }

    // right
    remainingDistance = expectedDistance;
    for (uint32_t i = index; i + 1 < junctions.size(); ++i)
    {
        GenomicRegion junctionRegion(
            junctions[i].getRefNameRight(),
            std::min(junctions[i].getPositionRight(), junctions[i+1].getPositionLeft()),
            std::max(junctions[i].getPositionRight(), junctions[i+1].getPositionLeft())
        );
        if (alignsWithinRegion(remainingDistance, junctionRegion, secondRecord, junctions[i].getPositionRight()))
            return true;
        if (remainingDistance <= 0)
            break;
    }
    if (remainingDistance > 0)
    {
        Junction lastJunction = junctions[junctions.size() - 1];
        GenomicRegion junctionRegion(
            lastJunction.getRefNameRight(),
            std::min(lastJunction.getPositionRight(), lastJunction.getPositionRight() + lastJunction.getDirectionRight() * 1000),
            std::max(lastJunction.getPositionRight(), lastJunction.getPositionRight() + lastJunction.getDirectionRight() * 1000)
        );
        if (alignsWithinRegion(remainingDistance, junctionRegion, secondRecord, lastJunction.getPositionRight()))
            return true;
    }
    return false;
}

bool ReadTemplate::alignsWithinRegion(int & remainingDistance, GenomicRegion & junctionRegion, BamRecord & secondRecord, int referencePosition)
{
    if (secondRecord.getAlignmentRegion().overlaps(junctionRegion))
    {
        int distanceFromJunction = std::min(
            std::abs(referencePosition - secondRecord.getStartPos()),
            std::abs(referencePosition - secondRecord.getEndPos())
        );
        if (distanceFromJunction <= remainingDistance)
            return true;
    } else {
        remainingDistance -= junctionRegion.getRegionEnd() - junctionRegion.getRegionStart();
    }
    return false;
}



void ReadTemplate::findSplitsBasedOnGaps(std::vector<Junction> & junctions)
{
    //for (auto & junction : junctions) {
    for (uint32_t i = 0; i < junctions.size(); ++i) {
        auto & junction = junctions[i];
        if (junction.getRefNameLeft() != junction.getRefNameRight())
            continue;
        if (junction.getDirectionLeft() == junction.getDirectionRight())
            continue;

        if (!this->splitFirst && this->records[this->primaryFirst].containsDeletion()) {
            if (readSpansJunction(this->records[this->primaryFirst], junction))
            {
                if (alignsWithinExpectedDistance(junctions, i, this->records[this->primaryFirst], this->records[this->primaryLast]))
                {
                    this->splittingJunctions.insert(junction.getID());
                    this->splitFirst = true;
                }
            }
        }
        
        if (!this->splitLast && this->records[this->primaryLast].containsDeletion())
        {
            if (readSpansJunction(this->records[this->primaryLast], junction))
            {
                if (alignsWithinExpectedDistance(junctions, i, this->records[this->primaryLast], this->records[this->primaryFirst]))
                {
                    this->splittingJunctions.insert(junction.getID());
                    this->splitLast = true;
                }
            }
        }
    }
    this->split = this->split || (this->splitFirst || this->splitLast);
}

bool ReadTemplate::readSpansJunction(BamRecord & record, Junction & junction)
{
    int requiredOverlap = 20;
    int deletionMargin = 5;

    if (record.getReferenceName() != junction.getRefNameLeft())
        return false;
    
    int deletionSize = std::abs(junction.getPositionRight() - junction.getPositionLeft());
    
    bool overlapLeft = false;
    bool overlapRight = false;
    if (record.getDeletionSize() >= deletionSize - deletionMargin && record.getDeletionSize() <= deletionSize + deletionMargin)
    {
        if (record.getStartPos() <= (junction.getPositionLeft() - requiredOverlap) && junction.getDirectionLeft() < 0)
            overlapLeft = true;
        if (record.getEndPos() >= (junction.getPositionRight() + requiredOverlap) && junction.getDirectionRight() > 0)
            overlapRight = true;
    }
    return (overlapLeft && overlapRight);
}

void ReadTemplate::findSpanningReads(std::vector<Breakpoint> & breakpoints)
{
    if (! this->isProperPair() || this->insertSize > 1000 || this->insertSize < 0)
        return;
    if (this->orientation == "FF" || this->orientation == "RR")
        return;
    
    if (this->records[this->primaryFirst].getReferenceName() != this->records[this->primaryLast].getReferenceName())
        return;
    
    int requiredOverlap = 20;
    
    for (Breakpoint & breakpoint : breakpoints)
    {
        if (this->records[this->primaryFirst].getReferenceName() != breakpoint.getReferenceName())
            continue;
        if (! this->splitFirst)
        {
            if ((breakpoint.getPosition() - this->records[this->primaryFirst].getStartPos()) >= requiredOverlap && (this->records[primaryFirst].getEndPos() - breakpoint.getPosition()) >= requiredOverlap)
            {
                this->spanningBreakpoints.insert(breakpoint.getID());
                this->spanning = true;
            }
        }
        if (! this->splitLast)
        {
            if (breakpoint.getPosition() - this->records[this->primaryLast].getStartPos() >= requiredOverlap && this->records[this->primaryLast].getEndPos() - breakpoint.getPosition() >= requiredOverlap)
            {
                this->spanningBreakpoints.insert(breakpoint.getID());
                this->spanning = true;
            }
        }
    }

    return;
}


void ReadTemplate::findBridgedBreakpoints(std::vector<Breakpoint> & breakpoints)
{
    if (! this->isProperPair() || this->insertSize > 1000 || this->insertSize < 0)
        return;
    if (this->orientation == "FF" || this->orientation == "RR")
        return;
    
    if (this->records[this->primaryFirst].getReferenceName() != this->records[this->primaryLast].getReferenceName())
        return;

    for (Breakpoint & breakpoint : breakpoints)
    {
        if (!this->records[this->primaryFirst].isReverse() && this->records[this->primaryLast].isReverse())
        {
            if (this->records[this->primaryFirst].getEndPos() <= breakpoint.getPosition() && this->records[this->primaryLast].getStartPos() >= breakpoint.getPosition())
                this->bridgedBreakpoints.insert(breakpoint.getID());
        }
        if (this->records[this->primaryFirst].isReverse() && !this->records[this->primaryLast].isReverse())
        {
            if (this->records[this->primaryFirst].getStartPos() >= breakpoint.getPosition() && this->records[this->primaryLast].getEndPos() <= breakpoint.getPosition())
                this->bridgedBreakpoints.insert(breakpoint.getID());
        }
    }
    return;
}


int ReadTemplate::softClippingLeftLength(BamRecord & r)
{
    return r.getClipLeft();
}

int ReadTemplate::softClippingRightLength(BamRecord & r)
{
    return r.getClipRight();
}


std::string ReadTemplate::getOrientation()
{
    return this->orientation;
}

int64_t ReadTemplate::getInsertSize()
{
    return this->insertSize;
}

bool ReadTemplate::containsSplitRead()
{
    return this->split;
}

bool ReadTemplate::alignsAcrossChromosomes()
{
    return this->interChromosome;
}

std::string ReadTemplate::getPrimaryRecordString()
{
    std::string tempString = "";
    if (this->primaryFirst == -1 || this->primaryLast == -1)
        return tempString;
    
    tempString.append(std::to_string(this->records[this->primaryFirst].getStartPos()));
    tempString.append("\t");
    tempString.append(std::to_string(this->records[this->primaryFirst].getEndPos()));
    tempString.append("\t");
    tempString.append(std::to_string(this->records[this->primaryFirst].isReverse()));
    tempString.append("\n");
    tempString.append(std::to_string(this->records[this->primaryLast].getStartPos()));
    tempString.append("\t");
    tempString.append(std::to_string(this->records[this->primaryLast].getEndPos()));
    tempString.append("\t");
    tempString.append(std::to_string(this->records[this->primaryLast].isReverse()));

    return tempString;
}

bool ReadTemplate::isProperPair()
{
    return (this->primaryFirst != -1 && this->primaryLast != -1);
}

std::vector<BamRecord> & ReadTemplate::getRecords()
{
    return this->records;
}

std::vector<BamRecord> ReadTemplate::getPrimaryRecords()
{
    std::vector<BamRecord> primaryRecords;
    if (this->primaryFirst != -1)
        primaryRecords.push_back(this->records[this->primaryFirst]);
    if (this->primaryLast != -1)
        primaryRecords.push_back(this->records[this->primaryLast]);
    return primaryRecords;
}

bool ReadTemplate::isInterestingReadPair(ReadPairFilter & filter)
{
    std::string rName1, rName2;

    std::vector<BamRecord> records = getPrimaryRecords();
    if (records.size() != 2)
        return false;

    return filter.fragmentPassesFilter(records[0], records[1]);
}

bool ReadTemplate::containsSpanningRead()
{
    return this->spanning;
}

float ReadTemplate::getTemplateWeight()
{
    return this->templateWeight;
}

void ReadTemplate::print()
{
    if (this->isProperPair()){
        std::cout << this->records[this->primaryFirst].getStartPos() << "\t" << this->records[this->primaryFirst].getEndPos() << "\t";
        std::cout << this->records[this->primaryFirst].getClipLeft() << "\t" << this->records[this->primaryFirst].getClipRight() << "\t";
        std::cout << this->records[this->primaryFirst].getSeqLength() << std::endl;
        std::cout << this->records[this->primaryLast].getStartPos() << "\t" << this->records[this->primaryLast].getEndPos() << "\t";
        std::cout << this->records[this->primaryLast].getClipLeft() << "\t" << this->records[this->primaryLast].getClipRight() << "\t";
        std::cout << this->records[this->primaryLast].getSeqLength() << std::endl;
        std::cout << std::endl;
    }
}

void ReadTemplate::setGCContent(float gcContent)
{
    this->templateGCContent = gcContent;
}

float ReadTemplate::getGCContent()
{
    return this->templateGCContent;
}

std::string ReadTemplate::getName()
{
	return this->templateName;
}

std::string ReadTemplate::getJunctionString()
{
    return this->junctionString;
}

std::string ReadTemplate::getBpSpanString()
{
    return this->bpSpanString;
}

std::string ReadTemplate::getBpBridgeString()
{
    return this->bpBridgeString;
}

std::vector<int> ReadTemplate::getMappingQualities()
{
    std::vector<int> mapQ(2, 0);
    mapQ[0] = this->records[this->primaryFirst].getMapQ();
    mapQ[1] = this->records[this->primaryLast].getMapQ();
    return mapQ;
}

void ReadTemplate::markSuspectedSplit()
{
	this->suspectedSplit = true;
}

bool ReadTemplate::containsSuspectedSplit()
{
	return this->suspectedSplit;
}