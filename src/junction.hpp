#ifndef JUNCTIONHEADER
#define JUNCTIONHEADER

#include <string>
#include <vector>
#include <seqan/bam_io.h>
#include "bamFileHandler.hpp"
#include "breakpoint.hpp"
#include "record.hpp"

class Junction {
    int id;
    std::string variantRefName;
    std::string rNameLeft; 
    std::string rNameRight;
    int directionLeft; 
    int directionRight;
    int xLeft;
    int xRight;

    int beginLeft, endLeft;
    int beginRight, endRight;

    Junction();

    public:
    Junction(int, std::string, std::string, std::string,int, int,int, int);
    Junction(int, std::string, std::string, int, int, bool);
    void adjustCoordinatesForBAM();
    void setRightSide(std::string, int, int);
    void setLeftSide(std::string, int, int);
    Breakpoint toBreakpoint();
    Breakpoint leftSideToBreakpoint(int);
    Breakpoint rightSideToBreakpoint(int);
    bool operator==(Junction);
    

    int getID();
    int getDirectionLeft();
    int getDirectionRight();
    int getPositionLeft();
    int getPositionRight();
    std::string getRefNameLeft();
    std::string getRefNameRight();
    std::string getVariantRefName();

    void setId(int);
    void setDirectionLeft(int);
    void setDirectionRight(int);
    void setPositionRight(int);
    void setPositionLeft(int);
    void setRefNameLeft(std::string);
    void setRefNameRight(std::string);
    void setVariantRefName(std::string);
    bool hasValidRightSide();
    bool hasValidLeftSide();

    bool qualifiesAsBreakpoint();
    std::vector<std::string> getAltStrings();

    void print();
};

#endif
