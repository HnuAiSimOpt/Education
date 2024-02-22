#ifndef Element_H_
#define Element_H_
#include "IntArray.h"
#include "DoubleMatrix.h"
#include "Section.h"

class Element
{
protected:
	IntArray nodes;
	int SectionId;
public:
	void computeStiffness(DoubleMatrix& answer);
	void getAssembleIndex(IntArray& loc);
};
#endif

