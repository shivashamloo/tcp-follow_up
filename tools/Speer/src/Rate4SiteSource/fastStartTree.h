// $Id: fastStartTree.h,v 1.1 2008/11/17 16:45:31 lanczyck Exp $

#ifndef ___FAST_START_TREE
#define ___FAST_START_TREE

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include <iostream>

using namespace std;



tree getBestMLTreeFromManyNJtrees(sequenceContainer & allTogether,
								stochasticProcess& sp,
								const int numOfNJtrees,
								const MDOUBLE tmpForStartingTreeSearch,
								const MDOUBLE epslionWeights,
								ostream& out);


#endif
