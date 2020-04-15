// 	$Id: fromCountTableComponentToDistance2USSRV.cpp,v 1.1 2008/11/17 16:45:31 lanczyck Exp $	

#include "fromCountTableComponentToDistance2USSRV.h"
#include "likeDist.h"
#include <cassert>

fromCountTableComponentToDistance2USSRV::fromCountTableComponentToDistance2USSRV(
		const countTableComponentGam& ctcBase,
		const countTableComponentHom& ctcSSRV,
		const ussrvModel &model,
		MDOUBLE toll,
		MDOUBLE brLenIntialGuess ) : _model(model), _ctcBase(ctcBase), _ctcSSRV(ctcSSRV) {
	_distance = brLenIntialGuess ;//0.03;
	_toll = toll;
}

void fromCountTableComponentToDistance2USSRV::computeDistance() {
	likeDist2USSRV likeDist1(_model,_toll); 
	MDOUBLE initGuess = _distance;
	_distance = likeDist1.giveDistance(_ctcBase,_ctcSSRV,_likeDistance,initGuess);
	assert(_distance>=0);
}
