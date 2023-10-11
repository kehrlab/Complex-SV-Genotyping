#include "genomicRegion.hpp"


GenomicRegion::GenomicRegion()
{
    setReferenceName("");
    setRegionStart(0);
    setRegionEnd(0);
    setReverse(false);
    this->containsSequence = false;
}

GenomicRegion::GenomicRegion(std::string rName, int begin, int end)
{
    setReferenceName(rName);
    setRegionStart(begin);
    setRegionEnd(end);
    setReverse(false);
    this->containsSequence = false;
}

GenomicRegion::GenomicRegion(std::string rName, int begin, int end, bool reverse)
{
    setReferenceName(rName);
    setRegionStart(begin);
    setRegionEnd(end);
    setReverse(reverse);
    this->containsSequence = false;
}

GenomicRegion::GenomicRegion(std::string regionString)
{
    std::vector<std::string> parts1;
    std::vector<std::string> parts2;
    boost::algorithm::split(parts1, regionString, boost::is_any_of(":"));
    boost::algorithm::split(parts2, parts1[1], boost::is_any_of("-"));
    setReferenceName(parts1[0]);
    setRegionStart(std::stoi(parts2[0]));
    setRegionEnd(std::stoi(parts2[1]));
    setReverse(false);
    this->containsSequence = false;
}

void GenomicRegion::setRegionStart(int begin)
{
    this->begin = begin;
}

void GenomicRegion::setRegionEnd(int end)
{
    this->end = end;
}

void GenomicRegion::setReferenceName(std::string referenceName)
{
    this->rName = referenceName;
}

void GenomicRegion::setReverse(bool reverse)
{
    this->reverse = reverse;
}

std::string GenomicRegion::getReferenceName()
{
    return this->rName;
}

int GenomicRegion::getRegionStart()
{
    return this->begin;
}

int GenomicRegion::getRegionEnd()
{
    return this->end;
}

bool GenomicRegion::overlaps(GenomicRegion other)
{
    if (this->getReferenceName() == other.getReferenceName())
    {
        if (other.getRegionStart() >= this->getRegionStart() && other.getRegionStart() <= this->getRegionEnd())
            return true;
        if (other.getRegionEnd() >= this->getRegionStart() && other.getRegionEnd() <= this->getRegionEnd())
            return true;
        if (other.getRegionStart() <= this->getRegionStart() && other.getRegionEnd() >= this->getRegionEnd())
            return true;
    }
    return false;
}

bool GenomicRegion::overlaps(GenomicRegion other, int requiredOverlap)
{
    if (this->getReferenceName() == other.getReferenceName())
    {
        if (other.getRegionStart() >= this->getRegionStart() && other.getRegionStart() <= this->getRegionEnd() && (other.getRegionStart() - this->getRegionStart() >= requiredOverlap || this->getRegionEnd() - other.getRegionStart() >= requiredOverlap))
            return true;
        if (other.getRegionEnd() >= this->getRegionStart() && other.getRegionEnd() <= this->getRegionEnd() && (other.getRegionEnd() - this->getRegionStart() >= requiredOverlap || this->getRegionEnd() - other.getRegionEnd() >= requiredOverlap))
            return true;
        if (other.getRegionStart() <= this->getRegionStart() && other.getRegionEnd() >= this->getRegionEnd())
            return true;
    }
    return false;
}

void GenomicRegion::mergeWithRegion(GenomicRegion other)
{
    this->begin = std::min(this->getRegionStart(), other.getRegionStart());
    this->end = std::max(this->getRegionEnd(), other.getRegionEnd());
}

void GenomicRegion::print()
{
    std::cout << getReferenceName() << "\t" << getRegionStart() << "\t" << getRegionEnd() << "\t" << (isReverse() ? "Reverse" : "Forward") << std::endl;
}

bool GenomicRegion::operator==(GenomicRegion rhs)
{
    return (
        rhs.getReferenceName() == this->getReferenceName() && 
        rhs.getRegionStart() == this->getRegionStart() && 
        rhs.getRegionEnd() == this->getRegionEnd()
    );
}

std::vector<GenomicRegion> GenomicRegion::joinRegions(std::vector<GenomicRegion> regions)
{
    for (auto it = regions.begin(); it != regions.end();)
    {
        bool merged = false;
        for (auto it1 = it+1; it1 != regions.end(); )
        {
            if (it->overlaps(*it1))
            {
                it->mergeWithRegion(*it1);
                it1 = regions.erase(it1);
                merged = true;
                break;
            } else {
                ++it1;
            }
        }
        if (!merged)
            ++it;
    }
    return regions;
}

void GenomicRegion::readSequence(SeqFileHandler & refFileHandler)
{
    this->sequence = refFileHandler.getRegionSequence(getReferenceName(), getRegionStart(), getRegionEnd());
    if (this->reverse)
        this->sequence = seqan::reverseComplementString(this->sequence);
    if (seqan::length(this->sequence) > 0)
        this->containsSequence = true;
}

void GenomicRegion::reverseRegion()
{
    if (this->reverse)
        this->reverse = false;
    else
        this->reverse = true;
    this->sequence = seqan::reverseComplementString(this->sequence);
}

seqan::Dna5String & GenomicRegion::getSequence()
{
    return this->sequence;
}

bool GenomicRegion::isReverse()
{
    return this->reverse;
}

void GenomicRegion::determineNStats(SeqFileHandler & refFileHandler)
{
    bool sequencePresent = false;
    if (!this->containsSequence)
        readSequence(refFileHandler);
    else
        sequencePresent = true;
    int nN = 0;
    int headRun = 0;
    int tailRun = 0;
    int longest = 0;
    int current = 0;
    seqan::Dna5String nString = 'N';
    bool begin = true;
    for (seqan::Dna5String s : this->sequence)
    {
        if (s == nString)
        {
            current += 1;
        } else {
            if (begin) {
                headRun = current;
                begin = false;
            }
            longest = std::max(longest, current);
            nN += current;
            current = 0;
        }
    }
    if (!sequencePresent)
        this->sequence = "";
    if (current > 0) {
        tailRun = current;
        nN += current;
        longest = std::max(longest, current);
        current = 0;
    }
    this->nFraction = (float) nN / seqan::length(this->sequence);
    this->longestN = longest;
    this->nStart = headRun;
    this->nEnd = tailRun;
}

bool GenomicRegion::isValidSample(SeqFileHandler & refFileHandler)
{
    determineNStats(refFileHandler);
    if (this->nFraction > 0.01)
        return false;
    if (this->longestN > 5)
        return false;
    if (this->nStart > 20 || this->nEnd > 20)
        return false;
    return true;
}

bool GenomicRegion::sequenceInMemory()
{
    return this->containsSequence;
}

void GenomicRegion::clearSequence()
{
    if (this->containsSequence)
    {
        this->sequence = "";
        this->containsSequence = false;
    }
}

std::string GenomicRegion::getRegionString()
{
    return (this->rName + ":" + std::to_string(this->begin) + "-" + std::to_string(this->end));
}