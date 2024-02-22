#pragma once
#include <vector>
#include <memory>
#include "Node.h"
#include "Element.h"

class LinearStatic
{
protected:
	DoubleMatrix Stiffness;
	DoubleArray force;
	DoubleArray disp;

	Domain* d;
public:
	LinearStatic(Domain* dd) : d(dd){}
	void assemble();
	
};
