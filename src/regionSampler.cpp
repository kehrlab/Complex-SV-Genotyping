#include "regionSampler.hpp"

RegionSampler::RegionSampler()
{}

RegionSampler::RegionSampler(std::vector<ContigInfo> contigInfo, std::vector<GenomicRegion> regions)
{
    this->gcBias = false;
    this->intervalSizes = 25000;
    this->totalIntervalLength = 2000000;
    this->contigInfo = contigInfo;
    
    sampleInsertSizeRegions(regions);
}

void RegionSampler::sampleInsertSizeRegions(std::vector<GenomicRegion> regions)
{
    determineCommonContigs();
    std::vector<GenomicRegion> wholeGenomeRegions = regions;
    if (wholeGenomeRegions.size() == 0 || !allRegionsValid(wholeGenomeRegions)) 
	    wholeGenomeRegions = createGenomicRegions();
    
    subsampleRegions(wholeGenomeRegions);
}

void RegionSampler::determineCommonContigs()
{
    for (unsigned i = 0; i < this->contigInfo[0].cNames.size(); ++i)
    {
        int minContigLength = this->contigInfo[0].cLengths[i];
        bool include = true;
        for (unsigned j = 0; j < this->contigInfo.size(); ++j)
        {
            bool present = false;
            for (unsigned k = 0; k < this->contigInfo[j].cNames.size(); ++k)
            {
                if (this->contigInfo[j].cNames[k] == this->contigInfo[0].cNames[i])
                {
                    present = true;
                    minContigLength = std::min(minContigLength, this->contigInfo[j].cLengths[k]);
                    break;
                }
            }
            if (!present)
            {
                include = false;
                break;
            }
        }
        if (this->contigInfo[0].cNames[i].size() > 5)
            include = false;
            
        if (include)
        {
            this->commonContigs.push_back(this->contigInfo[0].cNames[i]);
            this->minContigLengths.push_back(minContigLength);
        }
    }
}

bool RegionSampler::allRegionsValid(std::vector<GenomicRegion> regions)
{
    bool valid = true;
    for (GenomicRegion& region : regions)
    {
        bool validRegion = false;
        for (int i = 0; i < this->commonContigs.size(); ++i) 
	    if (region.getReferenceName() == this->commonContigs[i] && region.getRegionEnd() < this->minContigLengths[i])
                validRegion = true;
        valid = valid && validRegion;
    }
    return valid;
}

std::vector<GenomicRegion> RegionSampler::createGenomicRegions()
{
    std::vector<GenomicRegion> regions;
    for (unsigned j = 0; j < this->commonContigs.size(); ++j)
        regions.push_back(GenomicRegion(this->commonContigs[j], 0, this->minContigLengths[j]));    
    return regions;
}

void RegionSampler::subsampleRegions(std::vector<GenomicRegion> wholeGenomeRegions)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    int currentSize = 0;
    std::uniform_int_distribution<> regionDistrib(0, wholeGenomeRegions.size() - 1);
    int j = 0;
    while (currentSize < this->totalIntervalLength)
    {
        ++j;
	
        GenomicRegion r = wholeGenomeRegions[regionDistrib(gen)];
	if (r.getReferenceName() == "chrM" || r.getReferenceName() == "M" || r.getReferenceName() == "MT")
		continue;
	if (r.getRegionEnd() - r.getRegionStart() <= this->intervalSizes)
		continue;
	
        std::uniform_int_distribution<> distrib(0, r.getRegionEnd() - r.getRegionStart() - this->intervalSizes);
        int start = r.getRegionStart() + distrib(gen);	
        GenomicRegion sampledRegion(r.getReferenceName(), start, start + this->intervalSizes);

        // check N
        // if (this->gcBias)
        //     sampledRegion.readSequence(referenceFileHandler);
        // if (!sampledRegion.isValidSample(referenceFileHandler))
        //     continue;
        
        // add
        currentSize += this->intervalSizes;
        this->insertSizeRegions.push_back(sampledRegion);
    }
    this->insertSizeRegions = GenomicRegion::joinRegions(this->insertSizeRegions); 
}

std::vector<GenomicRegion> & RegionSampler::sampledInsertRegions()
{
    return this->insertSizeRegions;
}
