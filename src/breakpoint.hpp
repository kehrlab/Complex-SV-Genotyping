#ifndef BREAKPOINTHEADER
#define BREAKPOINTHEADER

#include <string>

class Breakpoint {
    int id;
    int x;
    std::string rName;
    
    public:
    Breakpoint(int, int, std::string);
    bool operator==(Breakpoint);
    int getPosition();
    int getID();
    std::string getReferenceName();
    void print();

    static struct {
        bool operator()(Breakpoint & a, Breakpoint & b) {
            return (a.getPosition() < b.getPosition());
        }
    } compareBreakpoints;
};

#endif