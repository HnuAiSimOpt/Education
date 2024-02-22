#pragma once
#include <vector>

class Node;
class Element;
class Domain
{
protected:
	std::vector<Node*> nodeList;
	std::vector<Element*> elementList;


public:
	int getNumberOfElementList() const { return elementList.size(); }
	std::vector<Element*>& getElementList()  { return elementList; }
	void initial();
};
