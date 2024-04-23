#ifndef CUSTOM_TYPE_HEADER
#define CUSTOM_TYPE_HEADER

#include <seqan/sequence.h> 
#include <seqan/bam_io.h>
#include <string>
#include <unordered_map>
#include <vector>

#include "genomicRegion.hpp"
#include "junction.hpp"
#include "breakpoint.hpp"

struct ContigInfo
{
    std::vector<std::string> cNames;
    std::vector<int32_t> cLengths;
    std::unordered_map<std::string, uint64_t> globalPositions;

    ContigInfo(){}

    ContigInfo(std::unordered_map<std::string, int> contigLengths)
    {
        std::vector<std::string>().swap(this->cNames);
        std::vector<int>().swap(this->cLengths);
        for (auto c : contigLengths)
        {
            this->cNames.push_back(c.first);
            this->cLengths.push_back(c.second);
        }
        sortNames();
        calculateGlobalContigPositions();
    }

    std::unordered_map<std::string, int32_t> getContigLengths()
    {
        std::unordered_map<std::string, int32_t> contigLengths;
        for (uint32_t i = 0; i < this->cNames.size(); ++i)
            contigLengths[this->cNames[i]] = this->cLengths[i];
        return contigLengths;
    }
    
    void sortNames() {
        std::vector<std::pair<std::string, int>> cInfo;
        for (uint32_t i = 0; i < this->cLengths.size(); ++i)
            cInfo.push_back(std::pair(this->cNames[i], this->cLengths[i]));
        std::sort(
            cInfo.begin(), 
            cInfo.end(), 
            [](std::pair<std::string, int32_t> p1, std::pair<std::string, int32_t> p2) {
                    return p1.second >= p2.second;
                }
            );
        for (uint32_t i = 0; i < cInfo.size(); ++i) {
            this->cNames[i] = cInfo[i].first;
            this->cLengths[i] = cInfo[i].second;
        }
    }

    void calculateGlobalContigPositions() {
        uint64_t startPos = 0;
        for (uint64_t i = 0; i < this->cNames.size(); ++i) {
            this->globalPositions[this->cNames[i]] = startPos;
            startPos += this->cLengths[i];
        }
    }
};

struct TemplatePosition
{
    int begin, end;
    std::string chr;
};

struct JunctionRegion
{
    std::vector<GenomicRegion> regions;
    std::vector<int32_t> junctionIndices;
    std::vector<int32_t> breakpointIndices;
    std::vector<Junction> junctions;
    std::vector<Breakpoint> breakpoints;
    int32_t length;
    std::string chromosome;

    void print()
    {
        std::cout << "-------------------------------------------" << std::endl;
        for (auto & r : this->regions)
            r.print();
        std::cout << "-------------- Junctions ------------------" << std::endl;
        for (auto & j : this->junctions)
            j.print();
        for (auto & j : this->junctionIndices)
            std::cout << j << "\t";
        std::cout << std::endl;
        std::cout << "-------------- Breakpoints ------------------" << std::endl;
        for (auto & b : this->breakpoints)
            b.print();
        for (auto & b : this->breakpointIndices)
            std::cout << b << "\t";
        std::cout << std::endl;
        std::cout << "Length: " << length << std::endl;
    }
};

struct BreakpointRegion
{
    GenomicRegion region;
    std::vector<Breakpoint> breakpoints;
};

struct VariantRegions
{
    std::vector<JunctionRegion> regions;
    std::vector<int32_t> distanceFromLast;
    std::vector<int32_t> distanceToNext;
};

struct SplitAlignmentInfo
{
    std::vector<std::unordered_set<int>> junctionIndices;
    std::vector<int> insertSize; 
};

#endif
