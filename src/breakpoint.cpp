#include "breakpoint.hpp"

Breakpoint::Breakpoint(int id, int x, std::string rName)
{
    this->id = id;
    this->x = x;
    this->rName = rName;
}

int Breakpoint::getID()
{
    return this->id;
}

int Breakpoint::getPosition()
{
    return this->x;
}

std::string Breakpoint::getReferenceName()
{
    return this->rName;
}

void Breakpoint::print()
{
    std::cout << this->id << " :\t" << this->rName << "\t" << this->x << std::endl; 
}

bool Breakpoint::operator==(Breakpoint rs)
{
    return (this->getPosition() == rs.getPosition() && this->getReferenceName() == rs.getReferenceName());
}