#include "Domain.h"
#include "Node.h"

void Domain::initial()
{
	Node* inode = new Node(1, this);
	nodeList.push_back(inode);
}
