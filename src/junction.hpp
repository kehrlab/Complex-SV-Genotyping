#ifndef JUNCTIONHEADER
#define JUNCTIONHEADER

#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>

#include "breakpoint.hpp"
#include "record.hpp"


class Junction {
    int id;
    std::string variantRefName;
    std::string rNameLeft; 
    std::string rNameRight;
    int directionLeft; 
    int directionRight;
    int32_t xLeft;
    int32_t xRight;

    int32_t beginLeft, endLeft;
    int32_t beginRight, endRight;

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
    int32_t getPositionLeft();
    int32_t getPositionRight();
    std::string getRefNameLeft();
    std::string getRefNameRight();
    std::string getVariantRefName();

    void setId(int);
    void setDirectionLeft(int);
    void setDirectionRight(int);
    void setPositionRight(int32_t);
    void setPositionLeft(int32_t);
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
