
#include <ncbi_pch.hpp>
#include "zToPTable.hpp"
#include "zToPTablePrivate.hpp"


USING_NCBI_SCOPE;

////////////////////////////////////////////////////
//
//   CZ_To_P_Table Class
//
////////////////////////////////////////////////////

const double CZ_To_P_Table::zInitial = zInitialHardcoded;
const double CZ_To_P_Table::zDelta = zDeltaHardcoded;
const unsigned int CZ_To_P_Table::pArraySize = sizeof(pOfZ)/sizeof(double);

//  Not using pArraySize in case order of declaration of statics is not predictable.
const valarray<double> CZ_To_P_Table::pArray(pOfZ, sizeof(pOfZ)/sizeof(double));  


double CZ_To_P_Table::GetP(double z)
{
    double p = -1.0;
    double zAbs = (z >= 0) ? z : -z;
    int upperBoundInt = -1;
    unsigned int lowerBound = 0, upperBound = 0;

    GetBoundingIndices(zAbs, lowerBound, upperBoundInt);
    upperBound = (unsigned int) upperBoundInt;

    //  Off end of the table...
    if (upperBoundInt < 0) {
        p = pArray[pArraySize - 1];

    //  ...exactly one of the table entries
    } else if (upperBound == lowerBound) {
        p = pArray[lowerBound];

    //  ... interpolate
    } else if (upperBound == lowerBound + 1) {
        p = InterpolateP(zAbs, lowerBound);

    } else {
        cerr << "Should never get here.... " << upperBound << "; " << lowerBound << endl;
    }

    //  Get the value for z < 0 if necessary.
    if (z != zAbs) {
        p = 1.0 - p;
    }

    return p;
}

//  Return true if inside the range of table; false otherwise.
bool CZ_To_P_Table::IsInRange(double z)
{
    double zAbs = (z >= 0) ? z : -z;
    return (zAbs <= GetMaxZ());
}

//  Depends on how many entries are in the table.
double CZ_To_P_Table::GetMaxZ()
{
    double maxZ = zInitial + zDelta*pArraySize;
    cerr << "GetMaxZ:  " << maxZ << endl;
    return maxZ;
}

double CZ_To_P_Table::GetP_MaxZ()  //  P(Z_max)
{
    return pArray[pArraySize - 1];
}

double CZ_To_P_Table::GetMinZ()    //  this is a negative number = -(GetMaxZ())
{
    return -GetMaxZ();
}

double CZ_To_P_Table::GetP_MinZ()  //  P(Z_min)
{
    return 1.0 - GetP_MaxZ();
}

//  Only works with z >= 0.
//  If bound is < 0, the value for z is off the corresponding end of the table.
//  Returns false (and both bounds < 0) if z < 0.
bool CZ_To_P_Table::GetBoundingIndices(double z, unsigned int& lowerBound, int& upperBound)
{
    static const double invZDelta = 1.0/zDelta;

    if (z < 0) {
        lowerBound = 0;
        upperBound = -1;
        return false;
    }

    unsigned int uiLower = (unsigned int) (z * invZDelta);
    lowerBound = (int) uiLower;
    upperBound = lowerBound + 1;

    if ((unsigned int)upperBound > pArraySize - 1) {
        lowerBound = (int) pArraySize - 1;
        upperBound = -1;
    }

    return true;
}

//  Only works with z >= 0.  
//  To interpolate for z < 0, use property P(z) + P(-z) = 1.
//  P(z) = P[lowerBound] + slope*(z - zDelta*lowerBound)
double CZ_To_P_Table::InterpolateP(double z, unsigned int lowerBound)
{

    static const double invZDelta = 1.0/zDelta;

    double slope, pInterpolated;
    double excessZ = (z - lowerBound * zDelta);

    if (z < 0 || lowerBound >= pArraySize) {
        pInterpolated = 1.0;
    } else {
        if (lowerBound == pArraySize - 1) {
            slope = -pArray[lowerBound] * invZDelta;
        } else {
            slope = (pArray[lowerBound + 1] - pArray[lowerBound]) * invZDelta;
        }
        pInterpolated = pArray[lowerBound] + slope*excessZ;  //-pArray[lowerBound] * invZDelta;
    }

/*
    cerr << "Interpolate z = " << z << "; lower bound " << lowerBound << "; excess z = " << excessZ << endl;
    cerr << "    pArray[lb]   = " << pArray[lowerBound] << endl;
    if (lowerBound < pArraySize - 1) {
        cerr << "    pArray[lb+1] = " << pArray[lowerBound+1] << endl;
    }
    cerr << "    result = " << pInterpolated << endl;
*/

    return pInterpolated;
}

