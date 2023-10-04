#include "readTemplate.hpp"
#include "bamFileHandler.hpp"
#include "junction.hpp"
#include "record.hpp"
#include <string>

ReadTemplate::ReadTemplate()
{}

ReadTemplate::ReadTemplate(std::vector<BamRecord> records)
{
    this->records = records;
    this->split = false;
    this->spanning = false;
    this->insertSize = 0;
    this->interChromosome = false;
    this->templateWeight = 0;
    this->templateGCContent = 0;
    this->templateName = records[0].getTemplateName();
    this->regionString = "";
    this->junctionString = "";
    this->chromosomeString = "";
    this->bpString = "";
    this->splitFirst = false;
    this->splitLast = false;
    this->suspectedSplit = false;

    determineInformativeRecords();
    determineFivePrimeEnds();
    determineOrientation();
    calculateInsertSize();
}

ReadTemplate::ReadTemplate(std::vector<BamRecord> records, std::vector<Junction> & junctions, std::vector<Breakpoint> & breakpoints)
{
    this->records = records;
    this->split = false;
    this->spanning = false;
    this->insertSize = 0;
    this->interChromosome = false;
    this->templateWeight = 0;
    this->templateGCContent = 0;
    this->templateName = records[0].getTemplateName();
    this->regionString = "";
    this->junctionString = "";
    this->bpString = "";
    this->chromosomeString = "";
    this->splitFirst = false;
    this->splitLast = false;
    this->suspectedSplit = false;
    determineInformativeRecords();
    determineFivePrimeEnds();
    determineOrientation();
    calculateInsertSize();
    findSplitReads(junctions);
    findSpanningReads(breakpoints);
}

void ReadTemplate::determineInformativeRecords()
{
    this->primaryFirst = -1;
    this->primaryLast = -1;
    for (int i = 0; i < this->records.size(); ++i)
    {
        if (records[i].isFirst())
            this->firstIndices.push_back(i);
        else if (records[i].isLast())
            this->lastIndices.push_back(i);
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

void ReadTemplate::determineFivePrimeEnds()
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
}

void ReadTemplate::calculateInsertSize()
{
    this->insertSize = 0;

    if (isProperPair())
    {
        if (this->records[this->primaryFirst].getReferenceName() != this->records[this->primaryLast].getReferenceName())
        {
            this->insertSize = 0;
            this->interChromosome = true;
            std::vector<std::string> cNamePair;
            cNamePair.push_back(this->records[this->primaryFirst].getReferenceName());
            cNamePair.push_back(this->records[this->primaryLast].getReferenceName());
            this->interChromosomeNames.push_back(cNamePair);
            return;
        }
    
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
            // if (this->insertSize < 0)
            //     this->insertSize--;
            // else
            //     this->insertSize++;
            return;
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

void ReadTemplate::determineGroup(std::vector<Junction> & junctions, std::vector<Breakpoint> & breakpoints, std::vector<GenomicRegion> & regions)
{
    findSplitReads(junctions);
    findSpanningReads(breakpoints);
    determineOverlappingRegions(regions);
}

void ReadTemplate::determineOverlappingRegions(std::vector<GenomicRegion> & regions)
{
    GenomicRegion fR = this->records[this->primaryFirst].getAlignmentRegion();
    GenomicRegion lR = this->records[this->primaryLast].getAlignmentRegion();

    for (int i = 0; i < regions.size(); ++i) {
        if (fR.overlaps(regions[i], 10) && !this->records[this->primaryFirst].isClipped())
            this->overlappingRegions.insert(i);
        if (lR.overlaps(regions[i], 10) && !this->records[this->primaryLast].isClipped())
            this->overlappingRegions.insert(i);
    }
}

void ReadTemplate::determineLocationStrings()
{
    if (this->overlappingRegions.size() > 0)
    {
        std::vector<int> regions;
        for (auto id : this->overlappingRegions)
        regions.push_back(id);
        std::sort(regions.begin(), regions.end());
        this->regionString = "";
        for (int i : regions)
            this->regionString = this->regionString + std::to_string(i);
    }
    
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
        this->bpString = "";
        for (int i : breakpoints)
            this->bpString = this->bpString + std::to_string(i);

    }
}

void ReadTemplate::findSplitReads(std::vector<Junction> & junctions)
{
    // both reads must properly align
    if (!isProperPair())
        return;
    findSplitsBasedOnClipping(junctions);
    findSplitsBasedOnGaps(junctions);
    return;
}

void ReadTemplate::findSplitsBasedOnClipping(std::vector<Junction> & junctions)
{
    BamRecord rFirst = this->records[this->primaryFirst];
    BamRecord rLast = this->records[this->primaryLast];

    // only works with at least one clipped record 
    if (!rFirst.isClipped() && !rLast.isClipped())
	    return;

    // get unique chromosome names
    std::unordered_set<std::string> chromosomes;
    for (Junction & junction : junctions)
        chromosomes.insert(junction.getVariantRefName());

    // search on each chromosome
    for (std::string cName : chromosomes)
        findClippedSplitsOnChromosome(cName, junctions);
}

void ReadTemplate::findClippedSplitsOnChromosome(std::string cName, std::vector<Junction> & junctions)
{
    BamRecord rFirst = this->records[this->primaryFirst];
    BamRecord rLast = this->records[this->primaryLast];

    std::vector<Junction> chromosomeJunctions;
    for (Junction & junction : junctions)
        if (junction.getVariantRefName() == cName)
            chromosomeJunctions.push_back(junction);

    int clippedLength = -1;
    int direction = 0;
    // bool interChromosomeTrue = this->interChromosome;
    for (unsigned i = 0; i < chromosomeJunctions.size(); ++i)
    {
        Junction junction = chromosomeJunctions[i];
        if (!this->splitFirst && isClippedAtJunction(clippedLength, direction, rFirst, junction))
        {
            if (alignsWithinExpectedDistance(chromosomeJunctions, i, rFirst, rLast))
            {
                this->splittingJunctions.insert(junction.getID());
                this->splitFirst = true;
                // check for interChromosome reads
                // if (junction.getRefNameLeft() != junction.getRefNameRight())
                // {
                //     if (interChromosomeTrue)
                //     {
                //         this->interChromosomeNames.clear();
                //         interChromosomeTrue = false;
                //     }
                //     std::vector<std::string> cNamePair;
                //     cNamePair.push_back(junction.getRefNameLeft());
                //     cNamePair.push_back(junction.getRefNameRight());
                //     this->interChromosomeNames.push_back(cNamePair);
                //     this->interChromosome = true;
                // }
            }
        }
        if (!this->splitLast && isClippedAtJunction(clippedLength, direction, rLast, junction))
        {
            if (alignsWithinExpectedDistance(chromosomeJunctions, i, rLast, rFirst))
            {
                this->splittingJunctions.insert(junction.getID());
                this->splitLast = true;
                // check for interChromosome reads
                // if (junction.getRefNameLeft() != junction.getRefNameRight())
                // {
                //     if (interChromosomeTrue)
                //     {
                //         this->interChromosomeNames.clear();
                //         interChromosomeTrue = false;
                //     }
                //     std::vector<std::string> cNamePair;
                //     cNamePair.push_back(junction.getRefNameLeft());
                //     cNamePair.push_back(junction.getRefNameRight());
                //     this->interChromosomeNames.push_back(cNamePair);
                //     this->interChromosome = true;
                // }
            }
        }
    }
    this->split = this->splitFirst || this->splitLast;
    return;
}

bool ReadTemplate::alignsWithinExpectedDistance(std::vector<Junction> & junctions, int index, BamRecord & clippedRecord, BamRecord & secondRecord)
{
    int expectedDistance = 1500;
    int remainingDistance = expectedDistance;

    // left 
    for (int i = index; i > 0; --i)
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
    for (int i = index; i < junctions.size() - 1; ++i)
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
    for (int i = 0; i < junctions.size(); ++i) {
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
    
    int deletionSize = std::abs(junction.getPositionRight() - junction.getPositionLeft()); // Problem!
    
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

bool ReadTemplate::isClippedAtJunction(int & clippedLength, int & direction, BamRecord & r, Junction & junction)
{
    std::string rName = r.getReferenceName();

    GenomicRegion regionLeft;
    GenomicRegion regionRight;

    int requiredOverlap = 20;
    int largeMargin = 20;
    int smallMargin = 10;

    if (junction.getDirectionLeft() < 0)
    {
        regionLeft = GenomicRegion(junction.getRefNameLeft(), junction.getPositionLeft() - largeMargin, junction.getPositionLeft() + smallMargin);
        int clipLen = softClippingRightLength(r);
        if (clipLen >= requiredOverlap)
        {
            GenomicRegion readEnd(rName, r.getEndPos(), r.getEndPos());
            if (regionLeft.overlaps(readEnd))
            {
                clippedLength = r.getSeqLength() - clipLen;
                direction = -1;
                return true;
            }
        }
    } 
    else 
    {
        regionLeft = GenomicRegion(junction.getRefNameLeft(), junction.getPositionLeft() - smallMargin, junction.getPositionLeft() + largeMargin);
        int clipLen = softClippingLeftLength(r);
        if (clipLen >= requiredOverlap)
        {
            GenomicRegion readEnd(rName, r.getStartPos(), r.getStartPos());
            if (regionLeft.overlaps(readEnd))
            {
                clippedLength = r.getSeqLength() - clipLen;
                direction = -1;
                return true;
            }
        }
    }


    if (junction.getDirectionRight() < 0)
    {
        regionRight = GenomicRegion(junction.getRefNameRight(), junction.getPositionRight() - largeMargin, junction.getPositionRight() + smallMargin);
        int clipLen = softClippingRightLength(r);
        if (clipLen >= requiredOverlap)
        {
            GenomicRegion readEnd(rName, r.getEndPos(), r.getEndPos());
            if (regionRight.overlaps(readEnd))
            {
                clippedLength = r.getSeqLength() - clipLen;
                direction = 1;
                return true;
            }
        }
    }
    else
    {
        regionRight = GenomicRegion(junction.getRefNameRight(), junction.getPositionRight() - smallMargin, junction.getPositionRight() + largeMargin);
        int clipLen = softClippingLeftLength(r);
        if (clipLen >= requiredOverlap)
        {
            GenomicRegion readEnd(rName, r.getStartPos(), r.getStartPos());
            if (regionRight.overlaps(readEnd))
            {
                clippedLength = r.getSeqLength() - clipLen;
                direction = 1;
                return true;
            }
        }
    }


    clippedLength = -1;
    direction = 0;
    return false;
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

int ReadTemplate::getInsertSize()
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

std::vector<std::vector<std::string>> ReadTemplate::getInterChromosomeNames()
{
    return this->interChromosomeNames;
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

std::string ReadTemplate::getRegionString()
{
    return this->regionString;
}

std::string ReadTemplate::getJunctionString()
{
    return this->junctionString;
}

std::string ReadTemplate::getBreakpointString()
{
    return this->bpString;
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

std::string ReadTemplate::getChromosomeString()
{
    this->chromosomeString = "";
    if (!this->isProperPair())
    {
        return "";
    }
    else
    {
        if (this->interChromosomeNames.size() > 0)
        {
            if (this->interChromosomeNames[0].size() == 2)
            {
                std::sort(this->interChromosomeNames[0].begin(), this->interChromosomeNames[0].end());
                for (auto & chrName : this->interChromosomeNames[0])
                    this->chromosomeString += chrName;
            }
        }
    }
    return this->chromosomeString;
}