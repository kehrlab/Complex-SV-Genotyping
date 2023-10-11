#include "vcfWriter.hpp"

VcfWriter::VcfWriter()
{}

VcfWriter::VcfWriter(std::string filename)
{
    setFileName(filename);
    openVcfFile();
}

VcfWriter::VcfWriter(std::string filename, std::vector<complexVariant> & variants, std::vector<std::string> fileNames, std::vector<std::vector<GenotypeResult>> & results)
{
    setFileName(filename);
    getSampleNames(fileNames);
    getContigNames(variants);
    createRefIDs();
    initHeader();
    createVcfRecords(variants, results);
}

void VcfWriter::setFileName(std::string filename)
{
    this->vcfFileName = filename;
}

void VcfWriter::getContigNames(std::vector<complexVariant> & variants)
{
    std::unordered_set<std::string> cNames;
    for (complexVariant cSV : variants)
    {
        std::unordered_set<std::string> vcNames = cSV.getContigNames();
        for (std::string cName : vcNames)
            cNames.insert(cName);
    }
    for (std::string cName : cNames)
        this->contigNames.push_back(cName);
}

void VcfWriter::getSampleNames(std::vector<std::string> fileNames)
{
    for (std::string fileName : fileNames)
    {
        std::string sampleName;
        // remove path information
        int splitIndex = fileName.rfind("/");
        if (splitIndex == fileName.npos)
            int splitIndex = fileName.rfind("\\");
        if (splitIndex != fileName.npos)
            sampleName = fileName.substr(splitIndex + 1);

        // remove ending
        splitIndex = sampleName.rfind(".");
        if (splitIndex != sampleName.npos)
            sampleName = sampleName.substr(0, splitIndex);
        
        this->sampleNames.push_back(sampleName);
    }
}

void VcfWriter::initHeader()
{
    std::time_t time = std::time(nullptr);
    std::string ctimeStr(ctime(&time));
    ctimeStr.erase(std::remove(ctimeStr.begin(), ctimeStr.end(), '\n'), ctimeStr.cend());

    seqan::appendValue(this->header, seqan::VcfHeaderRecord("fileformat", "VCFv4.3"));
    seqan::appendValue(this->header, seqan::VcfHeaderRecord("fileDate", ctimeStr));
    seqan::appendValue(this->header, seqan::VcfHeaderRecord("source", "cSVGenotyperV0.9"));
    seqan::appendValue(this->header, seqan::VcfHeaderRecord("INFO", "<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"));
    seqan::appendValue(this->header, seqan::VcfHeaderRecord("INFO", "<ID=EVENT,Number=A,Type=String,Description=\"Name of the corresponding complex Variant\">"));
    seqan::appendValue(this->header, seqan::VcfHeaderRecord("INFO", "<ID=MATEID,Number=A,Type=String,Description=\"ID of breakend mate\">"));
    seqan::appendValue(this->header, seqan::VcfHeaderRecord("FILTER", "<ID=q20,Description=\"Quality below 20\">"));
    seqan::appendValue(this->header, seqan::VcfHeaderRecord("ID", "<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    seqan::appendValue(this->header, seqan::VcfHeaderRecord("ID", "<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"));
}

void VcfWriter::createVcfRecords(std::vector<complexVariant> & variants, std::vector<std::vector<GenotypeResult>> & results)
{
    int id = 0;
    for (int i = 0; i < variants.size(); ++i)
    {
        complexVariant cSV = variants[i];
        std::vector<GenotypeResult> variantResults = results[i];
        std::vector<VcfInfo> variantRecords;
        std::vector<Junction> junctions = cSV.getAllJunctions();
        for (Junction junction : junctions)
        {
            std::vector<std::string> altStrings = junction.getAltStrings();

            std::vector<std::string> alleles;
            for (Allele allele : cSV.getAlleles())
            {
                for (Junction alleleJunction : allele.getNovelJunctions())
                {
                    if (junction == alleleJunction) {
                        alleles.push_back(allele.getName());
                        break;
                    }
                }
            }
            
            VcfInfo leftInfo(
                this->rIDs[junction.getRefNameLeft()],
                junction.getPositionLeft(),
                "bnd_" + std::to_string(id),
                "N",
                altStrings[0],
                -1,
                "BND",
                cSV.getName(),
                "bnd_" + std::to_string(id + 1),
                alleles
            );
            ++id;
            VcfInfo rightInfo(
                this->rIDs[junction.getRefNameRight()],
                junction.getPositionRight(),
                "bnd_" + std::to_string(id),
                "N",
                altStrings[1],
                -1,
                "BND",
                cSV.getName(),
                "bnd_" + std::to_string(id - 1),
                alleles
            );
            ++id;

            variantRecords.push_back(leftInfo);
            variantRecords.push_back(rightInfo);
        }

        // infer allele counts
        for (int j = 0; j < variantRecords.size(); ++j)
            for (auto result : variantResults)
                variantRecords[j].extractGenotype(result);

        // merge records
        for (auto it = variantRecords.begin(); it != variantRecords.end(); ++it)
        {
            for (auto it1 = it + 1; it1 != variantRecords.end(); )
            {
                if (it->isSameBreakend(*it1))
                {
                    it->mergeWith(*it1);
                    std::string oldID = it1->getRecordID();
                    std::string newID = it->getRecordID();
                    it1 = variantRecords.erase(it1);
                    for (auto it2 = variantRecords.begin(); it2 != variantRecords.end(); ++it2)
                    {
                        if (it2->isSameBreakend(*it))
                            continue;
                        it2->updateMateID(oldID, newID);
                    }
                } else {
                    ++it1;
                }
            }
        }


        for (VcfInfo recordInfo : variantRecords)
            this->vcfRecords.push_back(recordInfo.getVcfRecord());
    }
    return;
}

void VcfWriter::write()
{
    openVcfFile();
    setFormalities();
    writeHeader();
    writeVcfRecords();
    closeVcfFile();
}

void VcfWriter::openVcfFile()
{
    seqan::open(this->vcfFile, this->vcfFileName.c_str());
    if (!seqan::isOpen(this->vcfFile))
        throw std::runtime_error("Could not open VCF file.");
}

void VcfWriter::createRefIDs()
{
    int rID = 0;
    for (std::string cName : this->contigNames)
    {
        this->rIDs[cName] = rID;
        ++rID;
    }
}

void VcfWriter::setContigNames()
{
    for (std::string cName : this->contigNames) 
        seqan::appendValue(seqan::contigNames(seqan::context(this->vcfFile)), seqan::CharString(cName.c_str()));
}

void VcfWriter::setSampleNames()
{
    for (std::string sampleName : this->sampleNames)
        seqan::appendValue(seqan::sampleNames(seqan::context(this->vcfFile)), seqan::CharString(sampleName.c_str()));
}

void VcfWriter::setFormalities()
{
    setContigNames();
    setSampleNames();
}

void VcfWriter::writeHeader()
{
    seqan::writeHeader(this->vcfFile, this->header);
}

void VcfWriter::writeVcfRecords()
{
    for (seqan::VcfRecord record : this->vcfRecords)
        seqan::writeRecord(this->vcfFile, record);
}

void VcfWriter::closeVcfFile()
{
    if (seqan::isOpen(this->vcfFile))
        seqan::close(this->vcfFile);
}