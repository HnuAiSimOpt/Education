#ifndef mathfem_H_
#define mathfem_H_
#include "DoubleMatrix.h"
//返回：如果值小于0，返回-1，否则返回1
inline double sgn(double i)
{
	return (i < 0. ? -1. : 1.);
}

//最大值
inline int max_fem(const int a, const int b) { return a >= b ? a : b; }
inline double max_fem(const double a, const double b) { return a >= b ? a : b; }

//最小值
inline int min_fem(const int a, const int b) { return a <= b ? a : b; }
inline double min_fem(const double a, const double b) { return a <= b ? a : b; }

/*
 *求解三次方程(a*x^3+b*x^2+c*x+d=0)
 *@r1:第一根
 *@r2:第二根
 *@r3:第三根
 *@nroot:根的数目
 */
void cubic(double a, double b, double c, double d, double *r1, double *r2, double *r3, int *nroot);
/**
 * \brief 求解对称方阵（2*2 or 3*3）的特征值和特征向量采用Jacobi法，计算精确且速度较慢
 * \param vec 特征向量
 * \param val 特征值
 * \param src 待求的对称矩阵(src存在拷贝，是否采用这个？？)
 */
void eigOfSymmtricMatrix(DoubleMatrix& vec, DoubleArray& val, DoubleMatrix src);
#endif