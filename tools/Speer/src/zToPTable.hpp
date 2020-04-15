#ifndef ZTOP_TABLE__HPP
#define ZTOP_TABLE__HPP

#include <valarray>
#include <corelib/ncbistd.hpp>

USING_NCBI_SCOPE;

class CZ_To_P_Table {

private:

    static const double zInitial;
    static const double zDelta;
    static const unsigned int pArraySize;
    static const valarray<double> pArray;

public:
    CZ_To_P_Table() {};
    ~CZ_To_P_Table() {};

    //  Performs a linear interpolation of the two nearest Z values.
    //  If z is outside range of table, interpolate between P(Z_max) and 0,
    //  or P(Z_min) and 1.0, depending on the sign of z.
    static double GetP(double z);

    //  Return true if inside the range of table (i.e., |z| <= GetMaxZ()); false otherwise.
    static bool IsInRange(double z);

    //  Depends on how many entries are in the table.
    static double GetMaxZ();
    static double GetP_MaxZ();  //  P(Z_max)

    static double GetMinZ();    //  this is a negative number = -(GetMaxZ())
    static double GetP_MinZ();  //  P(Z_min)

private:

    //  Only works with z >= 0.
    //  A returned upperBound < 0 means the value for z is off the high end of the table.
    //  lowerBound can never be negative if z>= 0.
    //  Returns false (and both bounds < 0) if z < 0.
    static bool GetBoundingIndices(double z, unsigned int& lowerBound, int& upperBound);

    //  Only works with z >= 0 (returns 1.0 if z or lowerBound are out of range).  
    //  To interpolate for z < 0, use property P(z) + P(-z) = 1.
    static double InterpolateP(double z, unsigned int lowerBound);

};

#endif  // ZTOP_TABLE__HPP
