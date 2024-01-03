#ifndef mathfem_H_
#define mathfem_H_
#include "DoubleMatrix.h"
//���أ����ֵС��0������-1�����򷵻�1
inline double sgn(double i)
{
	return (i < 0. ? -1. : 1.);
}

//���ֵ
inline int max_fem(const int a, const int b) { return a >= b ? a : b; }
inline double max_fem(const double a, const double b) { return a >= b ? a : b; }

//��Сֵ
inline int min_fem(const int a, const int b) { return a <= b ? a : b; }
inline double min_fem(const double a, const double b) { return a <= b ? a : b; }

/*
 *������η���(a*x^3+b*x^2+c*x+d=0)
 *@r1:��һ��
 *@r2:�ڶ���
 *@r3:������
 *@nroot:������Ŀ
 */
void cubic(double a, double b, double c, double d, double *r1, double *r2, double *r3, int *nroot);
/**
 * \brief ���ԳƷ���2*2 or 3*3��������ֵ��������������Jacobi�������㾫ȷ���ٶȽ���
 * \param vec ��������
 * \param val ����ֵ
 * \param src ����ĶԳƾ���(src���ڿ������Ƿ�����������)
 */
void eigOfSymmtricMatrix(DoubleMatrix& vec, DoubleArray& val, DoubleMatrix src);
#endif