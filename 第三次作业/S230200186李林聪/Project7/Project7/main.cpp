//S230200186 ���ִ�
//���������ε�Ԫ
//���ƽ��Ӧ������
//�ڵ�1��5��10Ϊ�����̶��ڵ�
//�ڵ�15�ܵ����µ�������СΪ5�����0.025�����ɱ�0.3��
//�����������ļ�����
//�������ⲿ��Eigen
#include <iostream>
#include <Eigen/Dense>
#include<stdio.h>
#include"Calculate.h"
using namespace Eigen;
using namespace std;
int main()
{
	MatrixXf Node(15, 3);
	Node << 1, 0, 0,
		2, 0.1, 0,
		3, 0.2, 0,
		4, 0.3, 0,
		5, 0, 0.1,
		6, 0.1, 0.1,
		7, 0.2, 0.1,
		8, 0.3, 0.1,
		9, 0.4, 0.1,
		10, 0, 0.2,
		11, 0.1, 0.2,
		12, 0.2, 0.2,
		13, 0.3, 0.2,
		14, 0.4, 0.2,
		15, 0.5, 0.2;
	//��Ԫ
	MatrixXf Ele(16, 4);
	Ele << 1, 1, 2, 6,
		2, 1, 6, 5,
		3, 2, 3, 7,
		4, 2, 7, 6,
		5, 3, 4, 8,
		6, 3, 8, 7,
		7, 4, 9, 8,
		8, 5, 6, 11,
		9, 5, 11, 10,
		10, 6, 7, 12,
		11, 6, 12, 11,
		12, 7, 8, 13,
		13, 7, 13, 12,
		14, 8, 9, 14,
		15, 8, 14, 13,
		16, 9, 15, 14;
	Calculate calculate(Node,Ele,2.1e7,0.025,0.3);
	calculate.initial();
	calculate.F(29) = -5.0;
	calculate.Displacement(0) = 0; calculate.Displacement(1) = 0; calculate.Displacement(8) = 0;
	calculate.Displacement(9) = 0; calculate.Displacement(18) = 0; calculate.Displacement(19) = 0;
	calculate.Solve();
	calculate.printf();
	system("pause");
	return 0;
}
