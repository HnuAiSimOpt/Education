#include "LinearStatic.h"



void LinearStatic::assemble()
{
	int numberofElement = d->getNumberOfElementList();
	for (auto& ielem : d->getElementList())
	{
		DoubleMatrix stiff;
		ielem->computeStiffness(stiff);
		
	}
}
