#pragma once
#include "Domain.h"

class basicObj
{
public:
	int number;
	Domain* d;
	basicObj(int _n, Domain* ii) : number(_n), d(ii){}
	~basicObj(){ }
};
