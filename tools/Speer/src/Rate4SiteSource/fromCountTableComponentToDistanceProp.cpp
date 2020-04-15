// $Id: fromCountTableComponentToDistanceProp.cpp,v 1.1 2008/11/17 16:45:31 lanczyck Exp $

#include "fromCountTableComponentToDistanceProp.h"
#include "likeDistProp.h"

fromCountTableComponentToDistanceProp::fromCountTableComponentToDistanceProp(
		const vector<countTableComponentGam>& ctc,
		const vector<stochasticProcess> &sp,
		const MDOUBLE toll,
		const MDOUBLE brLenIntialGuess ) : _sp(sp), _ctc(ctc) {
	_distance =brLenIntialGuess;
	_toll = toll;
}

void fromCountTableComponentToDistanceProp::computeDistance() {
	likeDistProp likeDist1(alphabetSize(),_sp,_toll);
	_distance = likeDist1.giveDistance(_ctc,_likeDistance);
}
