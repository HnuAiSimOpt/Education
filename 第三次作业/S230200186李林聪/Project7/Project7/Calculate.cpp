#include"Calculate.h"
Calculate::Calculate(MatrixXf node, MatrixXf element, float e, float m_t, float m_mu)
{
	this->Node = node;
	this->Element = element;
	this->E = e;
	this->t = m_t;
	this->mu = m_mu;
	this->num_Element = Element.rows();
	this->num_node = Node.rows();
	length = 2 * Node.rows();
	this->K = MatrixXf(num_node * 2, num_node * 2);
	K.setZero(num_node * 2, num_node * 2);
	F= VectorXf(length);
	F.setZero();
	Displacement= VectorXf(length);
	Displacement.setOnes();
}
void Calculate::initial()
{
	MatrixXf k(6, 6);
	for (int i = 0; i < num_Element; i++)
	{  
		k = Stiffness(Ele_Pos(i));
		Setserial(Element.block(i, 1, 1, 3));
		K = Assemble(K, k);
	}
}
MatrixXf Calculate::Ele_Pos(int i)
{
	MatrixXf temp(1, 3), b(1, 3), c(1, 3), d(1, 3);
	temp = Element.block(i, 1, 1, 3);             
	b = Node.block(temp(0, 0) - 1, 1, 1, 2);  
	c = Node.block(temp(0, 1) - 1, 1, 1, 2);  
	d = Node.block(temp(0, 2) - 1, 1, 1, 2);  
	MatrixXf node(3, 2);
	node << b, c, d;
	return node;
}
MatrixXf Calculate::Stiffness(MatrixXf coord)
{
	float A = (coord(0, 0) * (coord(1, 1) - coord(2, 1)) +
		coord(1, 0) * (coord(2, 1) - coord(0, 1)) + coord(2, 0) * (coord(0, 1) - coord(1, 1))) / 2;
	MatrixXf B = GetB(coord, A);
	MatrixXf D(3, 3);
	D << 1, mu, 0,
		mu, 1, 0,
		0, 0, (1 - mu) / 2;
	D = (E / (1 - mu * mu)) * D;
	MatrixXf k(6, 6);
	k = t * A * B.transpose() * D * B;
	return k;
}
MatrixXf Calculate::GetB(MatrixXf coord, float area)
{
	MatrixXf B(3, 6);
	B << coord(1, 1) - coord(2, 1), 0, coord(2, 1) - coord(0, 1), 0, coord(0, 1) - coord(1, 1), 0,
		0, coord(2, 0) - coord(1, 0), 0, coord(0, 0) - coord(2, 0), 0, coord(1, 0) - coord(0, 0),
		coord(2, 0) - coord(1, 0), coord(1, 1) - coord(2, 1), coord(0, 0) - coord(2, 0),
		coord(2, 1) - coord(0, 1), coord(1, 0) - coord(0, 0), coord(0, 1) - coord(1, 1);
	B = B / (2 * area);
	return B;
}
MatrixXf Calculate::Assemble(MatrixXf K, MatrixXf k)
{
	
	for (int n1 = 0; n1 < 6; n1++)
	{
		for (int n2 = 0; n2 < 6; n2++)
		{
			K(serial[n1], serial[n2]) = K(serial[n1], serial[n2]) + k(n1, n2);
		}
	}
	serial.clear();
	return K;
}
void Calculate::Setserial(MatrixXf num)
{
	serial.push_back(2 * num(0, 0) - 2);
	serial.push_back(2 * num(0, 0) - 1);
	serial.push_back(2 * num(0, 1) - 2);
	serial.push_back(2 * num(0, 1) - 1);
	serial.push_back(2 * num(0, 2) - 2);
	serial.push_back(2 * num(0, 2) - 1);
}
void Calculate::Solve()
{
	for (int j = 0; j < length; j++)
	{
		if (Displacement(j) == 0)
		{
			K.block(j, 0, 1, K.cols()).setZero();    
			K.block(0, j, K.rows(), 1).setZero();  
			K(j, j) = 1;
			F(j) = 0;    
		}
	}
	Displacement = K.lu().solve(F);         
	Force = VectorXf(length);
	Force = K * Displacement;
}
void Calculate::printf()
{
	cout << "rows of node:" << Node.rows() << endl;
	cout << "整体刚度矩阵 K：" << endl << K.block(Node.rows(), Node.rows(), Node.rows(), Node.rows()) << endl;
	cout << "节点位移：" << endl << Displacement << endl;
	cout << "支反力：\n" << Force << endl;
}
