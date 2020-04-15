// $Id: treeUtil.h,v 1.1 2008/11/17 16:45:33 lanczyck Exp $

#ifndef ___TREE_UTIL
#define ___TREE_UTIL
#include "definitions.h"
#include "tree.h"

vector<tree> getStartingTreeVecFromFile(string fileName);

tree starTree(const vector<string>& names);

void getStartingTreeVecFromFile(string fileName,
											vector<tree>& vecT,
											vector<char>& constraintsOfT0);


bool sameTreeTolopogy(tree t1, tree t2);

bool cutTreeToTwo(tree bigTree,
			  const string& nameOfNodeToCut,
			  tree &small1,
			  tree &small2);

tree::nodeP makeNodeBetweenTwoNodes(	tree& et,
										tree::nodeP nodePTR1,
										tree::nodeP nodePTR2,
										const string &interName);

void cutTreeToTwoSpecial(const tree& source,
						tree::nodeP intermediateNode,
						tree &resultT1PTR,
						tree &resultT2PTR);

vector<string> getSequencesNames(const tree& t);

#endif
