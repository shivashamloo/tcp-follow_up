// $Id: distanceTable.h,v 1.1 2008/11/17 16:45:30 lanczyck Exp $

#ifndef ___DISTANCE_TABLE
#define ___DISTANCE_TABLE

#include "definitions.h"
#include "distanceMethod.h"
#include "sequenceContainer.h"

void giveDistanceTable(const distanceMethod* dis,
					   const sequenceContainer& sc,
					   VVdouble& res,
					   vector<string>& names,
					   const vector<MDOUBLE> * weights = NULL);


#endif
