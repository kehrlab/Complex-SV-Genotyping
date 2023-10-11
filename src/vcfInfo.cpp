#include "vcfInfo.hpp"

VcfInfo::VcfInfo()
{
    setRefID(-1);
    setBeginPos(0);
    setRecordID(".");
    setRefString("N");
    this->alt = ".";
    setQuality(-1);
    setSVType("BND");
    setFormat();
    setFilter("PASS");
}

VcfInfo::VcfInfo(int rID, int beginPos, std::string id, std::string ref, std::string alt, int qual, std::string svtype, std::string event, std::string mateid, std::vector<std::string> alleles)
{
    setRefID(rID);
    setBeginPos(beginPos);
    setRecordID(id);
    setRefString(ref);
    setAltString(alt, event, mateid);
    setQuality(qual);
    setSVType(svtype);
    setAssociatedAlleles(alleles);
    setFormat();
    setFilter("PASS");
}

void VcfInfo::setRefID(int rID)
{
    this->rID = rID;
}

void VcfInfo::setBeginPos(int beginPos)
{
    this->beginPos = beginPos;
}

void VcfInfo::setRecordID(std::string id)
{
    this->id = id;
}

void VcfInfo::setRefString(std::string ref)
{
    this->ref = ref;
}

void VcfInfo::setAltString(std::string alt, std::string event, std::string mateid)
{
    this->altAlleles.push_back(alt);
    this->events.push_back(event);
    this->mateids.push_back(mateid);
}

void VcfInfo::setQuality(int qual)
{
    this->qual = qual;
}

void VcfInfo::setSVType(std::string svtype)
{
    this->svtype = svtype;
}

void VcfInfo::setAssociatedAlleles(std::vector<std::string> associatedAlleles)
{
    this->associatedAlleles = associatedAlleles;
}

void VcfInfo::setFormat()
{
    this->format = "GT:GQ";
}

void VcfInfo::setFilter(std::string filter)
{
    this->filter = filter;
}

void VcfInfo::createAltString()
{
    if (this->altAlleles.size() == 0)
        return;
    
    this->alt = "";
    for (int i = 0; i < this->altAlleles.size() - 1; ++i)
        this->alt = this->alt + this->altAlleles[i] + ",";
    this->alt = this->alt + this->altAlleles[this->altAlleles.size() - 1];
}

void VcfInfo::extractGenotype(GenotypeResult & result)
{
    std::string gt = result.getCalledGenotype();
    int q = result.getQuality();

    int splitIndex = gt.find("/");
    std::string allele1 = gt.substr(0, splitIndex);
    std::string allele2 = gt.substr(splitIndex + 1);

    int gtCount = 0;
    for (std::string allele : this->associatedAlleles) {
        if (allele == allele1)
            gtCount++;
        if (allele == allele2)
            gtCount++;
    }
    std::vector<int> alleleCounts;
    std::vector<float> alleleQu;
    alleleCounts.push_back(gtCount);
    alleleQu.push_back(q);
    this->genotypes.push_back(alleleCounts);
    this->gtQualities.push_back(alleleQu);
}

bool VcfInfo::isSameBreakend(VcfInfo recordInfo)
{
    bool sameLocation = (this->rID == recordInfo.getRefID() && this->beginPos == recordInfo.getBeginPos() && this->ref == recordInfo.getRefString() && this->events[0] == recordInfo.getEvents()[0]);
    bool overlap = false;

    if (this->altAlleles[0].substr(0,1) == "N" && recordInfo.getAltStrings()[0].substr(0, 1) == "N")
        overlap = true;
    if (this->altAlleles[0].substr(this->altAlleles[0].size() - 1) == "N" && recordInfo.getAltStrings()[0].substr(recordInfo.getAltStrings()[0].size() - 1) == "N")
        overlap = true;
    return (sameLocation && overlap);
}

void VcfInfo::mergeWith(VcfInfo recordInfo)
{
    std::string old_id = recordInfo.getRecordID();
    std::vector<std::string> altStrings = recordInfo.getAltStrings();
    std::vector<std::string> events = recordInfo.getEvents();
    std::vector<std::string> mateIDs = recordInfo.getMateIDs();

    std::vector<std::vector<int>> genotypes = recordInfo.getGenotypes();
    std::vector<std::vector<float>> gtQs = recordInfo.getGtQualities();

    for (int i = 0; i < altStrings.size(); ++i)
    {
        this->altAlleles.push_back(altStrings[i]);
        this->events.push_back(events[i]);
        this->mateids.push_back(mateIDs[i]);
        
        for (int j = 0; j < genotypes.size(); ++j)
        {
            this->genotypes[j].push_back(genotypes[j][i]);
            this->gtQualities[j].push_back(gtQs[j][i]);
        }
    }
}

void VcfInfo::updateMateID(std::string oldID, std::string newID)
{
    for (int i = 0; i < this->mateids.size(); ++i)
        if (this->mateids[i] == oldID)
            this->mateids[i] = newID;
}

int VcfInfo::getRefID()
{
    return this->rID;
}

int VcfInfo::getBeginPos()
{
    return this->beginPos;
}

std::string VcfInfo::getRecordID()
{
    return this->id;
}

std::string VcfInfo::getRefString()
{
    return this->ref;
}

std::vector<std::string> VcfInfo::getAltStrings()
{
    return this->altAlleles;
}

std::vector<std::string> VcfInfo::getEvents()
{
    return this->events;
}

std::vector<std::string> VcfInfo::getMateIDs()
{
    return this->mateids;
}

std::string VcfInfo::getFilter()
{
    return this->filter;
}

std::string VcfInfo::getSVType()
{
    return this->svtype;
}

std::vector<std::vector<int>> VcfInfo::getGenotypes()
{
    return this->genotypes;
}

std::vector<std::vector<float>> VcfInfo::getGtQualities()
{
    return this->gtQualities;
}

seqan::VcfRecord VcfInfo::getVcfRecord()
{
    // gather information
    createAltString();
    createInfoString();
    inferGenotypeFromAlleleCounts();

    // set record info
    seqan::VcfRecord record;
    record.rID = this->rID;
    record.beginPos = this->beginPos;
    record.alt = this->alt.c_str();
    record.ref = this->ref.c_str();
    record.qual = record.MISSING_QUAL();
    record.filter = this->filter.c_str();
    record.format = this->format.c_str();
    record.id = this->id.c_str();

    for (int i = 0; i < this->genotype.size(); ++i) 
        appendValue(record.genotypeInfos, (this->genotype[i] + ":" + std::to_string(this->genotypeQuality[i])).c_str());
    record.info = this->info.c_str();

    return record;
}

void VcfInfo::createInfoString()
{
    this->info = "SVTYPE=" + this->svtype + ";EVENT=";
    for (int i = 0; i < this->events.size(); ++i)
    {
        this->info = this->info + this->events[i];
        if (i < this->events.size() - 1)
            this->info = this->info + ",";
    }
    this->info = this->info + ";MATEID=";
    for (int i = 0; i < this->mateids.size(); ++i)
    {
        this->info = this->info + this->mateids[i];
        if (i < this->mateids.size() - 1)
            this->info = this->info + ",";
    }
}

void VcfInfo::inferGenotypeFromAlleleCounts()
{
    std::vector<std::vector<int>> alleleIndices;
    std::vector<int> alleleCounts;

    for (int j = 0; j < this->genotypes.size(); ++j)
    {
        std::vector<int> sampleAlleleIndices;
        for (int i = 0; i < this->genotypes[j].size(); ++i) {
            int alleleCount = this->genotypes[j][i];
            alleleCounts.push_back(alleleCount);
            if (alleleCount > 0) {
                sampleAlleleIndices.push_back(i);
            }
        }
        alleleIndices.push_back(sampleAlleleIndices);
    }
    
    for (int j = 0; j < alleleIndices.size(); ++j)
    {
        float minQuality = this->gtQualities[j][0];
        for (int i = 1; i < this->gtQualities[j].size(); ++i)
            minQuality = std::min(minQuality, this->gtQualities[j][i]);
        if (alleleIndices[j].size() == 0)
        {
            this->genotype.push_back("0/0");
        } else if (alleleIndices[j].size() == 1)
        {
            if (this->genotypes[j][alleleIndices[j][0]] == 2)
                this->genotype.push_back(std::to_string(alleleIndices[j][0]+1) + "/" + std::to_string(alleleIndices[j][0]+1));
            else if (this->genotypes[j][alleleIndices[j][0]] == 1)
                this->genotype.push_back("0/" + std::to_string(alleleIndices[j][0]+1));
            else
                std::cout << "Too many alleles" << std::endl;
        } else if (alleleIndices[j].size() == 2)
        {
            if (this->genotypes[j][alleleIndices[j][0]] == 1 && this->genotypes[j][alleleIndices[j][1]] == 1)
                this->genotype.push_back(std::to_string(alleleIndices[j][0]+1) + "/" + std::to_string(alleleIndices[j][1]+1));
            else
                std::cout << "Too many alleles" << std::endl;
        }
        this->genotypeQuality.push_back(minQuality);
    }
}