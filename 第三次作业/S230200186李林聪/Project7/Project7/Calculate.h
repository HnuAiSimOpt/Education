#pragma once
#include<iostream>
#include <Eigen/Dense>
#include<vector>
using namespace Eigen;
using namespace std;
class Calculate
{
    public:
	MatrixXf Node;
	MatrixXf Element;
	MatrixXf K;
	VectorXf F;
	VectorXf Displacement;
	VectorXf Force;
	int num_node;
	int num_Element;
	int length;
	float E;
	float t;
	float mu;
	Calculate(MatrixXf node, MatrixXf element,float m_e,float m_t,float m_mu);
	void initial();
	MatrixXf Ele_Pos(int i);
	MatrixXf Stiffness(MatrixXf ele_xy);
	MatrixXf Assemble(MatrixXf K, MatrixXf k);
	void Solve();
	void printf();
	MatrixXf GetB(MatrixXf coord, float area);
	std::vector<int>serial;
	void Setserial(MatrixXf num);
};

