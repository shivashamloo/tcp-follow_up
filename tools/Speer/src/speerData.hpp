#ifndef SPEER_DATA__HPP
#define SPEER_DATA__HPP

#include <string>
#include "speerTypes.hpp"

const unsigned int nAAReal = 20;
const double invNAAReal = 1.0/nAAReal;
const unsigned int nAA = 24;
const double invNAA = 1.0/nAA;
enum EAminoAcid {
    eAA_A = 0,
    eAA_C,
    eAA_D,
    eAA_E,
    eAA_F,
    eAA_G,
    eAA_H,
    eAA_I,
    eAA_K,
    eAA_L,
    eAA_M,
    eAA_N,
    eAA_P,
    eAA_Q,
    eAA_R,
    eAA_S,
    eAA_T,
    eAA_V,
    eAA_W,
    eAA_Y,
    eAA_IndeterminateB,  //  for "residue" 'B' (D or N)
    eAA_IndeterminateZ,  //  for "residue" 'Z' (E or Q)
    eAA_IndeterminateX,  //  for "residue" 'X' (any real AA)
    eAA_Gap              //  for the gap character '-' 
};
/*
struct EAAGroupCompare : public binary_function< EAAGroup, EAAGroup, bool > {
    bool operator()(const EAAGroup& lhs, const EAAGroup& rhs) const {
        return (int) lhs < (int) rhs;
    }
};
*/
const std::string::value_type AminoAcidChar[nAA] = { 'A',
                                      'C',
                                      'D',
                                      'E',
                                      'F',
                                      'G',
                                      'H',
                                      'I',
                                      'K',
                                      'L',
                                      'M',
                                      'N',
                                      'P',
                                      'Q',
                                      'R',
                                      'S',
                                      'T',
                                      'V',
                                      'W',
                                      'Y',
                                      'B',
                                      'Z',
                                      'X',
                                      '-'};
const std::string AA_B = "DN";
const std::string AA_Z = "EQ";
const std::string AA_X = "ACDEFGHIKLMNPQRSTVWY"; 

//  Background amino acid frequency map.
extern TAADataMap refFreqMap;

//  Map between a property name and the AA->property value map.
extern TPropertyMap speerPropertyMap;

//  Precomputed (unweighted) distances between all residue pairs:
//  sqrt{sum_m[xi_m - xj_m]}
//  where x_m is the property value associated with the i/j amino acid.
extern TPairDistanceMap speerPairDistanceMap;

//  The above maps are initialized with default values in this method.
//  It is possible that the main program may allow the user to input
//  additional properties and/or alternative values for pre-defined properties.
void SpeerInitializeDefaultGlobalData();


#endif  //  SPEER_DATA__HPP
