#include "junction.hpp"

Junction::Junction(
    int id, 
    std::string variantRefName, std::string rNameLeft, std::string rNameRight, 
    int xLeft, int xRight, 
    int directionLeft, int directionRight
) {
    setId(id);
    setVariantRefName(variantRefName);
    setLeftSide(rNameLeft, xLeft, directionLeft);
    setRightSide(rNameRight, xRight, directionRight);
}

Junction::Junction(
    int id,
    std::string variantRefName, 
    std::string rName, 
    int x, 
    int direction, 
    bool left
)
{
    this->rNameLeft = "";
    this->rNameRight = "";
    this->xLeft = -1;
    this->xRight = -1;
    this->directionLeft = 0;
    this->directionRight = 0;

    setId(id);
    setVariantRefName(variantRefName);
    if (left)
        setLeftSide(rName, x, direction);
    else
        setRightSide(rName, x, direction);
}

void Junction::setId(int id)
{
    this->id = id;
}

void Junction::setVariantRefName(std::string refName)
{
    this->variantRefName = refName;
}

void Junction::setRightSide(
    std::string rNameRight,
    int xRight, int directionRight
)
{
    setRefNameRight(rNameRight);
    setPositionRight(xRight);
    setDirectionRight(directionRight);
}

void Junction::setLeftSide(
    std::string rNameLeft,
    int xLeft, int directionLeft
)
{
    setRefNameLeft(rNameLeft);
    setPositionLeft(xLeft);
    setDirectionLeft(directionLeft);
}

Breakpoint Junction::toBreakpoint()
{
    if (!qualifiesAsBreakpoint())
    {
        throw std::runtime_error("ERROR: Conversion of Junction to Breakpoint failed: Junctions does not qualify as breakpoint.");
    }
    return Breakpoint(this->id, this->xLeft, this->rNameLeft);
}

bool Junction::qualifiesAsBreakpoint()
{
    return (
        this->xLeft == this->xRight && 
        this->rNameLeft == this->rNameRight &&
        this->directionLeft == -1 &&
        this->directionRight == 1
    );
}

Breakpoint Junction::leftSideToBreakpoint(int id)
{
    return Breakpoint(id, this->xLeft, this->rNameLeft);
}

Breakpoint Junction::rightSideToBreakpoint(int id)
{
    return Breakpoint(id, this->xRight, this->rNameRight);
}

void Junction::setDirectionLeft(int directionLeft)
{
    this->directionLeft = directionLeft;
}

void Junction::setDirectionRight(int directionRight)
{
    this->directionRight = directionRight;
}

void Junction::setPositionRight(int positionRight)
{
    this->xRight = positionRight;
}

void Junction::setPositionLeft(int positionLeft)
{
    this->xLeft = positionLeft;
}

void Junction::setRefNameLeft(std::string rNameLeft)
{
    this->rNameLeft = rNameLeft;
}

void Junction::setRefNameRight(std::string rNameRight)
{
    this->rNameRight = rNameRight;
}

void Junction::adjustCoordinatesForBAM()
{
    -- this->xLeft;
    -- this->xRight;
}

std::string Junction::getVariantRefName()
{
    return this->variantRefName;
}

int Junction::getDirectionLeft()
{
    return this->directionLeft;
}

int Junction::getDirectionRight()
{
    return this->directionRight;
}

int Junction::getPositionLeft()
{
    return this->xLeft;
}

int Junction::getPositionRight()
{
    return this->xRight;
}

std::string Junction::getRefNameLeft()
{
    return this->rNameLeft;
}

std::string Junction::getRefNameRight()
{
    return this->rNameRight;
}

int Junction::getID()
{
    return this->id;
}

void Junction::print()
{
    std::cout << this->id << " (" << this->variantRefName << ")" << " : " << std::endl;
    std::cout << "\tLeft :\t" << this->rNameLeft << "\t" << this->xLeft << "\t" << this->directionLeft << std::endl;
    std::cout << "\tRight :\t" << this->rNameRight << "\t" << this->xRight << "\t" << this->directionRight << std::endl;  
}

bool Junction::operator==(Junction rhs)
{
    if (this->rNameLeft != rhs.getRefNameLeft())
        return false;
    if (this->rNameRight != rhs.getRefNameRight())
        return false;
    if (this->xLeft != rhs.getPositionLeft())
        return false;
    if (this->xRight != rhs.getPositionRight())
        return false;
    if (this->directionLeft != rhs.getDirectionLeft())
        return false;
    if (this->directionRight != rhs.getDirectionRight())
        return false;
    if (this->variantRefName != rhs.getVariantRefName())
        return false;
    return true;
}

std::vector<std::string> Junction::getAltStrings()
{
    std::vector<std::string> altStrings(2, "");
    std::string pos = getRefNameRight() + ":" + std::to_string(getPositionRight());
    if (getDirectionRight() < 0)
        altStrings[0] = "N]" + pos + "]"; 
    else if (getDirectionRight() > 0)
        altStrings[0] = "N[" + pos + "[";
    pos = getRefNameLeft() +  ":" + std::to_string(getPositionLeft());
    if (getDirectionLeft() < 0)
        altStrings[1] = "]" + pos + "]N";
    else if (getDirectionLeft() > 0)
        altStrings[1] = "[" + pos + "[N";

    return altStrings;
}

bool Junction::hasValidLeftSide()
{
    return (this->rNameLeft != "" && this->xLeft >= 0 && this->directionLeft != 0);
}

bool Junction::hasValidRightSide()
{
    return (this->rNameRight != "" && this->xRight >= 0 && this->directionRight != 0);
}