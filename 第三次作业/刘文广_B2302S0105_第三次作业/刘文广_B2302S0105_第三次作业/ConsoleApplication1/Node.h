#ifndef Node_H_
#define Node_H_
#include "DoubleArray.h"
#include "basicObj.h"

class Node : public basicObj
{
protected:
	DoubleArray coords;
public:
	Node(int _n, Domain* ii);

	
};
#endif
