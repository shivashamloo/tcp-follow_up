// $Id: evaluateCharacterFreq.h,v 1.1 2008/11/17 16:45:31 lanczyck Exp $

#ifndef __Evaluate_Character_Freq_h
#define __Evaluate_Character_Freq_h

#include <iostream>
using namespace std;

#include "sequenceContainer.h"
#include "definitions.h"


vector<MDOUBLE> evaluateCharacterFreq(const sequenceContainer & sc);
VVdouble evaluateCharacterFreqOneForEachGene(const vector<sequenceContainer> & scVec);
vector<MDOUBLE> evaluateCharacterFreqBasedOnManyGenes(const vector<sequenceContainer> & scVec);

void changeCountsToFreqs(vector<MDOUBLE>& charFreq);
void makeSureNoZeroFreqs(vector<MDOUBLE> & charFreq);

#endif
