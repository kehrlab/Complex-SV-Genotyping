#ifndef CUSTOM_TYPE_HEADER
#define CUSTOM_TYPE_HEADER
#include <seqan/sequence.h> 
#include <seqan/bam_io.h>
#include <string>
#include <vector>

struct ContigInfo
{
    std::vector<std::string> cNames;
    std::vector<int> cLengths;
};

struct TemplatePosition
{
    int begin, end;
    std::string chr;
};

#endif
