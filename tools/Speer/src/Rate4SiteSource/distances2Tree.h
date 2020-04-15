// $Id: distances2Tree.h,v 1.1 2008/11/17 16:45:30 lanczyck Exp $

#ifndef ___DISTANCES2TREE
#define ___DISTANCES2TREE

#include "definitions.h"
#include "tree.h"
#include <string>
using namespace std;

class distances2Tree {
public:
  virtual ~distances2Tree() {}
  virtual  distances2Tree* clone() const =0;
  virtual tree computeTree(VVdouble distances, const vector<string>& names, const tree * const constriantTree = NULL) = 0;
};

#endif
