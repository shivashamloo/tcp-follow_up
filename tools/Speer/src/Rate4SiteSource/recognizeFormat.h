// $Id: recognizeFormat.h,v 1.1 2008/11/17 16:45:32 lanczyck Exp $

#ifndef ___RECOGNIZE_FORMAT
#define ___RECOGNIZE_FORMAT

#include "sequenceContainer.h"

class recognizeFormat{
public:
	static sequenceContainer read(istream &infile, const alphabet* alph);
	static void write(ostream &out, const sequenceContainer& sd);
	//readUnAligned: the input sequences do not need to be aligned (not all sequences are the same length).
	static sequenceContainer readUnAligned(istream &infile, const alphabet* alph);
};

#endif



