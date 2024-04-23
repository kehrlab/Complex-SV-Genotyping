#ifndef VARIANTPARSERHEADER
#define VARIANTPARSERHEADER

#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>

#include "json.hpp"

#include "junction.hpp"


typedef std::vector<std::vector<Junction>> variantData;

class variantParser {
    std::vector<variantData> variants;
    std::string filename;
    nlohmann::json jsonObject;
    int id;
    std::string refName;
    std::vector<std::vector<std::string>> alleleNames;
    std::vector<std::string> variantNames;

    void readStructure();
    void readJSONFile();
    void parseJSONObject();
    bool structureIsValid(std::vector<Junction>);

    variantParser();
    Junction parseJunction(nlohmann::detail::iteration_proxy_value<nlohmann::detail::iter_impl<nlohmann::basic_json<>>> &);

    public:
    variantParser(std::string filename);
    std::vector<variantData> getVariantJunctions();
    std::vector<std::vector<std::string>> getAlleleNames();
    std::vector<std::string> getVariantNames();
};

#endif
