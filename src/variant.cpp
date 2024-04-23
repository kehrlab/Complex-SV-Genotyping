#include "variant.hpp"

complexVariant::complexVariant()
{
    this->filterMargin = 500;
}

complexVariant::complexVariant(std::string variantName, std::vector<std::string> alleleNames, variantData variantJunctions)
{
    this->variantName = variantName;
    this->variantFileName = "";
    this->filterMargin = 500;

    for (uint32_t i = 0; i < variantJunctions.size(); ++i)
        this->variantAlleles.push_back(Allele(variantJunctions[i], alleleNames[i]));
    this->nAlleles = alleleNames.size();

    poolBreakpoints();
    poolJunctions();

    this->variantAlleles.push_back(Allele(this->allBreakpoints));
    calculateVariantRegions();
}

complexVariant::complexVariant(std::string variantName, std::vector<std::string> alleleNames, variantData variantJunctions, std::string variantFileName)
{
    this->variantName = variantName;
    this->variantFileName = variantFileName;
    this->filterMargin = 500;

    for (uint32_t i = 0; i < variantJunctions.size(); ++i)
        this->variantAlleles.push_back(Allele(variantJunctions[i], alleleNames[i]));
    this->nAlleles = alleleNames.size();

    poolBreakpoints();
    poolJunctions();

    this->variantAlleles.push_back(Allele(this->allBreakpoints));
    calculateVariantRegions();
}

complexVariant::complexVariant(std::vector<Junction> novelJunctions)
{
    this->variantName = "Unnamed Variant";
    this->variantFileName = "Unknown";
    this->nAlleles = 1;
    this->filterMargin = 500;
    this->variantAlleles.push_back(Allele(novelJunctions, "VAR"));
    poolBreakpoints();
    poolJunctions();
    calculateVariantRegions();
}


void complexVariant::print()
{
    std::cout << this->variantName << std::endl;
    std::cout << "----------------------" << std::endl;
    std::cout << "----------------------" << std::endl;
    for (auto it : this->variantAlleles)
        it.print();
    std::cout << "----------------------" << std::endl;
}

void complexVariant::poolBreakpoints()
{
    this->allBreakpoints.erase(this->allBreakpoints.begin(), this->allBreakpoints.end());
    for (auto allele : this->variantAlleles)
    {
        std::vector<Breakpoint> breakpoints = allele.getReferenceBreakpoints();
        for (auto bp : breakpoints)
        {
            bool present = false;
            for (auto it : this->allBreakpoints)
            {
                if (bp == it)
                {
                    present = true;
                    break;
                }
            }
            if (!present)
                this->allBreakpoints.push_back(bp);
        }
    }
}

void complexVariant::poolJunctions()
{
    this->allJunctions.erase(this->allJunctions.begin(), this->allJunctions.end());
    for (auto allele : this->variantAlleles)
    {
        std::vector<Junction> junctions = allele.getNovelJunctions();
        for (auto junction : junctions)
        {
            bool present = false;
            for (auto it : this->allJunctions)
            {
                if (junction == it)
                {
                    present = true;
                    break;
                }
            }
            if (!present)
                this->allJunctions.push_back(junction);
        }
    }
}

std::vector<GenomicRegion> complexVariant::calculateAssociatedRegions(LibraryDistribution & dist)
{
    double searchDistance = calculateSearchDistance(dist.getInsertMean(), dist.getInsertSD());
    std::vector<GenomicRegion> associatedRegions = createRegionsFromBreakpoints(searchDistance);
    associatedRegions = GenomicRegion::joinRegions(associatedRegions);
    return associatedRegions;
}

int complexVariant::calculateSearchDistance(double insertMean, double insertSD)
{
    return (int) (this->filterMargin + (insertMean + 3 * insertSD));
}

std::vector<GenomicRegion> complexVariant::createRegionsFromBreakpoints(int searchDistance)
{
    std::vector<GenomicRegion> associatedRegions;
    for (auto it : this->allBreakpoints)
        associatedRegions.push_back(
            GenomicRegion(it.getReferenceName(), std::max(it.getPosition() - searchDistance, 0), it.getPosition() + searchDistance)
        );
    return associatedRegions;
}

std::string complexVariant::getName()
{
    return this->variantName;
}

std::vector<Junction> & complexVariant::getAllJunctions()
{
    return this->allJunctions;
}

std::vector<Breakpoint> & complexVariant::getAllBreakpoints()
{
    return this->allBreakpoints;
}

std::vector<Allele> & complexVariant::getAlleles()
{
    return this->variantAlleles;
}

std::unordered_set<std::string> complexVariant::getContigNames()
{
    std::unordered_set<std::string> cNames;
    for (Junction junction : getAllJunctions())
    {
        cNames.insert(junction.getVariantRefName());
        cNames.insert(junction.getRefNameLeft());
        cNames.insert(junction.getRefNameRight());
    }
    return cNames;
}

void complexVariant::calculateVariantRegions()
{
    std::unordered_set<std::string> chromosomes;
    for (Breakpoint & bp : this->allBreakpoints)
        chromosomes.insert(bp.getReferenceName());
    
    for (std::string cName : chromosomes)
    {
        std::vector<Breakpoint> chrBreakpoints;

        for (Breakpoint & bp : this->allBreakpoints)
            if (bp.getReferenceName() == cName)
                chrBreakpoints.push_back(bp);
        
        if (chrBreakpoints.size() == 0)
            continue;

        this->variantRegions.push_back(GenomicRegion(cName, 0, chrBreakpoints[0].getPosition() - 1));
        for (uint32_t i = 0; i + 1 < chrBreakpoints.size(); ++i)
            this->variantRegions.push_back(GenomicRegion(cName, chrBreakpoints[i].getPosition(), chrBreakpoints[i+1].getPosition() - 1));
        this->variantRegions.push_back(GenomicRegion(cName, chrBreakpoints[chrBreakpoints.size() - 1].getPosition(), chrBreakpoints[chrBreakpoints.size() - 1].getPosition() + 100000));
    }
}

std::vector<GenomicRegion> & complexVariant::getVariantRegions()
{
    return this->variantRegions;
}

void complexVariant::createAlleleMaps(int filterMargin, std::unordered_map<std::string, int> contigLengths)
{
    for (auto & allele : this->variantAlleles)
        allele.createChromosomeMaps(this->allBreakpoints, filterMargin, contigLengths);
}

std::string complexVariant::getVariantFileName()
{
    return this->variantFileName;
}

void complexVariant::setFilterMargin(int filterMargin)
{
    this->filterMargin = filterMargin;
}
