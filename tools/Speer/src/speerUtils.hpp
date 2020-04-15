/*  $Id: speerUtils.hpp,v 1.2 2008/10/01 16:53:58 lanczyck Exp $  */

#ifndef SPEER_UTILS__HPP
#define SPEER_UTILS__HPP

#include <map>
#include <string>
#include <vector>
#include <corelib/ncbistd.hpp>
#include <corelib/ncbi_limits.hpp>

#include "speerTypes.hpp"
#include "speerAlignment.hpp"

USING_NCBI_SCOPE;
//USING_SCOPE(objects);
//USING_SCOPE(cd_utils);

//  Generate sequence weights a la Henicoff/Henicoff for each row of 
//  the input alignment.  The output map is indexed by the zero-based
//  row number.  If no columns passed the maxGapFraction, the map
//  'weights' will be returned empty.
//  Gaps, B, Z, X are treated as normal residues for simplicity in computing sequence weights.
void ComputeSequenceWeights(const CSpeerAlignment& alignment, TWeightMap& weights, double maxGapFraction = 0.20);

#endif // SPEER_UTILS__HPP
