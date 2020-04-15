// $Id: njConstrain.h,v 1.1 2008/11/17 16:45:32 lanczyck Exp $

#ifndef ___NJ_CONSTRAINT
#define ___NJ_CONSTRAINT

#include <map>


#include "sequenceContainer.h"
#include "tree.h"
using namespace std;

class njConstraint {
public:
  njConstraint(const tree& starttree, const tree& constraintTree);
  bool isCompatible(const tree::nodeP& n1, const tree::nodeP& n2, const bool verbose=false) const;
  void join(const tree::nodeP& n1, const tree::nodeP& n2, const tree::nodeP& newFather);
  void output(ostream &out) const;
  
private:
  tree _cTree;			// constriant tree
  map<tree::nodeP,tree::nodeP> _interTreeMap;


};

ostream &operator<<(ostream &out, const njConstraint &c);

#endif // ___NJ_CONSTRAINT
