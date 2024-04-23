#include "variantParser.hpp"

variantParser::variantParser(std::string filename)
{
    this->filename = filename;
    readStructure();
}


void variantParser::readStructure()
{
    readJSONFile();
    parseJSONObject();
}

void variantParser::readJSONFile()
{
    std::ifstream inStream(this->filename.c_str());
    try
    {
        inStream >> this->jsonObject;
    }
    catch (nlohmann::detail::parse_error const & e)
    {
        throw std::runtime_error("ERROR: Could not parse variant JSON file");
    }
}

void variantParser::parseJSONObject()
{
    for (auto& variant : this->jsonObject.items())
    {
        this->id = 0;
	bool isValid = true;
        variantData variantJunctions;
        std::vector<std::string> variantAlleleNames;
	std::string varName = variant.key();

        for (auto& allele : variant.value().items())
        {
            variantAlleleNames.push_back(allele.key());
            std::vector<Junction> tempJunctions;
            for (auto& chromosome : allele.value().items())
            {
                this->refName = chromosome.key();
                for (auto& junction : chromosome.value().items())
                {
                    ++ this->id;
                    tempJunctions.push_back(parseJunction(junction));
                }
            }
	    
	    // check whether the junctions for this allele form valid
            // chromosome structures
            if (!structureIsValid(tempJunctions)) {
                isValid = false;
                break;
            }

            variantJunctions.push_back(tempJunctions);
        }

        // do not add variant if its structure is not valid
        if (!isValid) {
            std::cerr << "Variant " << varName << " has invalid structure and will be ignored." << std::endl;
            continue;
        }	

	this->variantNames.push_back(varName);
        this->alleleNames.push_back(variantAlleleNames);
        this->variants.push_back(variantJunctions);
    }
}

bool variantParser::structureIsValid(std::vector<Junction> junctions)
{
	for (unsigned i = 0; i < junctions.size() - 1; ++i)
	{
		if (junctions[i].getVariantRefName() != junctions[i+1].getVariantRefName())
				continue;
		if (junctions[i].getRefNameRight() == junctions[i+1].getRefNameLeft())
		{
			if (junctions[i].getPositionRight() > junctions[i+1].getPositionLeft() && (junctions[i].getDirectionRight() > 0 || junctions[i+1].getDirectionLeft() < 0))
			{
				std::cerr << "ERROR: Junction positions indicate inverted segment but directions indicate normal orientation." << std::endl;
				return false;
			} else if (junctions[i].getPositionRight() <= junctions[i+1].getPositionLeft() && (junctions[i].getDirectionRight() < 0 || junctions[i+1].getDirectionLeft() > 0))
			{
				std::cerr << "ERROR: Junction positions indicate normal orientation but directions indicate inverted segment." << std::endl;
				return false;
			} else if (junctions[i].getPositionRight() == junctions[i+1].getPositionLeft())
			{
				std::cerr << "ERROR: Right breakpoint position of one junction matches left breakpoint position of next junction." << std::endl;
				return false;
			}
		}
	}
	return true;
}

Junction variantParser::parseJunction(nlohmann::detail::iteration_proxy_value<nlohmann::detail::iter_impl<nlohmann::basic_json<>>> & junction)
{
    std::string rNameLeft = "";
    std::string rNameRight = "";
    int xLeft = -1;
    int xRight = -1;
    int directionLeft = 0;
    int directionRight = 0;

    if (junction.value().find("rNameLeft") != junction.value().end())
        rNameLeft = junction.value()["rNameLeft"];
    if (junction.value().find("rNameRight") != junction.value().end())
        rNameRight = junction.value()["rNameRight"];
    if (junction.value().find("xLeft") != junction.value().end())
        xLeft = junction.value()["xLeft"];
    if (junction.value().find("xRight") != junction.value().end())
        xRight = junction.value()["xRight"];
    if (junction.value().find("directionLeft") != junction.value().end())
        directionLeft = (junction.value()["directionLeft"] == "right") ? 1 : -1;
    if (junction.value().find("directionRight") != junction.value().end())
        directionRight = junction.value()["directionRight"] == "right" ? 1 : -1;

    Junction tempJunction(
        this->id,
        this->refName,
        rNameLeft, rNameRight,
        xLeft - 1, xRight - 1,
        directionLeft, directionRight
    );
    return tempJunction;
}

std::vector<variantData> variantParser::getVariantJunctions()
{
    return this->variants;
}

std::vector<std::vector<std::string>> variantParser::getAlleleNames()
{
    return this->alleleNames;
}

std::vector<std::string> variantParser::getVariantNames()
{
    return this->variantNames;
}
