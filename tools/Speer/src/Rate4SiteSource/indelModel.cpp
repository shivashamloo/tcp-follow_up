// 	$Id: indelModel.cpp,v 1.1 2008/11/17 16:45:31 lanczyck Exp $	
#include "indelModel.h"


void indelModel::setFreqX(const MDOUBLE freq_x)
{
	_freq[0] =freq_x ;
	_alpha = 1/(2*_freq[0]*_freq[1]) ;
}

void indelModel::setFreqG(const MDOUBLE freq_g)
{
	_freq[0] =freq_g ;
	_alpha = 1/(2*_freq[0]*_freq[1]) ;
}
