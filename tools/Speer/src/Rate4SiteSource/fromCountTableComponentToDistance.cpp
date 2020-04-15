// $Id: fromCountTableComponentToDistance.cpp,v 1.1 2008/11/17 16:45:31 lanczyck Exp $

#include "fromCountTableComponentToDistance.h"
#include "likeDist.h"
#include <cassert>

fromCountTableComponentToDistance::fromCountTableComponentToDistance(
		const countTableComponentGam& ctc,
		const stochasticProcess &sp,
		const MDOUBLE toll,
		const MDOUBLE brLenIntialGuess ) : _sp(sp), _ctc(ctc) {
	_distance =brLenIntialGuess ;//0.03;
	_toll = toll;
}

void fromCountTableComponentToDistance::computeDistance() {
	likeDist likeDist1(_sp,_toll);
	MDOUBLE initGuess = _distance;
	_distance = likeDist1.giveDistance(_ctc,_likeDistance,initGuess);
	assert(_distance>=0);
}
