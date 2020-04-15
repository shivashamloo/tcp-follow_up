/*  $Id: speerTypes.hpp,v 1.3 2008/10/01 16:53:58 lanczyck Exp $  */

#ifndef SPEER_TYPES__HPP
#define SPEER_TYPES__HPP

#include <map>
#include <string>

using namespace std;

typedef map<unsigned int, double> TWeightMap;
typedef TWeightMap::iterator TWeightMapIt;
typedef TWeightMap::const_iterator TWeightMapCit;

//  These maps are written in terms of string::value_type to emphasize that 
//  we need to specify the character-type assumed in the definition of string.
typedef map< string::value_type, double > TAADataMap;
typedef TAADataMap::iterator TAADataMapIt;
typedef TAADataMap::const_iterator TAADataMapCit;
typedef TAADataMap::value_type TAADataVT;

//  This type is used to hold an arbitrary collection of maps between
//  amino acid residue type and quantitative properties of the amino acids.
typedef map<string, TAADataMap> TPropertyMap;
typedef TPropertyMap::iterator TPropertyMapIt;
typedef TPropertyMap::const_iterator TPropertyMapCit;
typedef TPropertyMap::value_type TPropertyMapVT;

typedef map<string::value_type, unsigned int> TCountsMap;
typedef TCountsMap::iterator TCountsMapIt;
typedef TCountsMap::const_iterator TCountsMapCit;

enum ESpeerCaseSensitivity {
    eSpeerCaseSensitive = 0,    //  honor differences between upper/lower case in alignment
    eSpeerConvertToUpper,       //  modify input strings to use all upper-case characters
    eSpeerConvertToLower        //  modify input strings to use all lower-case characters
};



typedef double TSpeerScore;
typedef pair< string::value_type, string::value_type > TResiduePair;


typedef map<unsigned int, TSpeerScore> TColumnScoreMap;
typedef TColumnScoreMap::iterator TColumnScoreMapIt;
typedef TColumnScoreMap::const_iterator TColumnScoreMapCit;

typedef map<TResiduePair, TSpeerScore> TPairDistanceMap;
typedef TPairDistanceMap::iterator TPairDistanceMapIt;
typedef TPairDistanceMap::const_iterator TPairDistanceMapCit;
typedef TPairDistanceMap::value_type TPairDistanceMapVT;


#endif // SPEER_TYPES__HPP
