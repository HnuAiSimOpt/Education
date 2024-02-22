#ifndef DOUBLEMATRIX_H
#define DOUBLEMATRIX_H
#include <vector>
#include <cmath>
#include "InfoOut.h"
using namespace std;
class DoubleArray;
class IntArray;
class DoubleMatrix
{
protected:
	//行数
	int nRows;
	//列数
	int nColumns;
	//值
	std::vector<double>values;
public:
	DoubleMatrix(int _nRows = 0, int _nColumns = 0) : nRows(_nRows), nColumns(_nColumns), values(_nRows*_nColumns)
	{
	}
	//拷贝构造函数
	DoubleMatrix(const DoubleMatrix& mat) : nRows(mat.nRows), nColumns(mat.nColumns), values(mat.values) {}
	DoubleMatrix(DoubleMatrix&& mat) : nRows(mat.nRows), nColumns(mat.nColumns), values(std::move(mat.values)){}
	DoubleMatrix& operator=(DoubleMatrix&& mat)
	{
		nRows = mat.nRows;
		nColumns = mat.nColumns;
		values = std::move(mat.values);
		return *this;
	}
	//赋值函数
	DoubleMatrix& operator=(const DoubleMatrix& mat)
	{
		nRows = mat.nRows;
		nColumns = mat.nColumns;
		values = mat.values;
		return *this;
	}
	~DoubleMatrix();
	//基本功能
	void SetValues(const int i, const int j, const double value);
	void SetnRows(const int i) { nRows = i; }
	void SetnColumns(const int j) { nColumns = j; }
	double GetValues(const int i, const int j) const { return (*this)(i, j); }
	int GetRow() const { return nRows; }
	int GetColumn() const { return nColumns; }
	void clear() { values.clear(); nRows = 0; nColumns = 0; }
	void resize(int rows, int columns) { nRows = rows; nColumns = columns; values.resize(nRows*nColumns); }
	//列表初始化
	DoubleMatrix(std::initializer_list<std::initializer_list<double>>src);
	DoubleMatrix& operator=(std::initializer_list<std::initializer_list<double>>src);
	//判断矩阵是否存在有
	inline bool isFinite() const
	{
		for (double val : values)
		{
			if (!std::isfinite(val))
			{
				return false;
			}
		}
		return true;
	}
	//判断是否为方阵
	inline bool isSquare() const
	{
		return nRows == nColumns;
	}
	//判断是否为空
	inline bool isEmpty() const
	{
		return nRows == 0 || nColumns == 0;
	}
	//得到i行j列元素
	//可修改的左值
	inline double& operator() (const int& i, const int& j) { return values[nColumns * (i - 1) + j - 1]; }
	inline double operator() (const int& i, const int& j) const { return values[nColumns * (i - 1) + j - 1]; }
	double& at(int i, int j) { return values[nColumns * i + j]; }
	double at(int i, int j) const { return values[nColumns * i + j]; }
	//求矩阵行列式的值
	double ComputeDeterminant() const;
	//矩阵运算
	//注意不是相同维度也能相加
	void Add(const DoubleMatrix& src, double s = 1);
	void minus(const DoubleMatrix& src);
	void Scale(const double& cof);
	DoubleMatrix operator+(const DoubleMatrix &src) const;
	DoubleMatrix operator-(const DoubleMatrix &src) const;
	DoubleMatrix operator*(const DoubleMatrix &src) const;
	DoubleMatrix operator*(const double& x) const;
	DoubleMatrix Inverse() const;
	void beInverseOf(const DoubleMatrix& src);
	void beInverseAndTransposeOf(const DoubleMatrix& src);
	DoubleMatrix Tanspose() const;
	void beTanspose() = delete;
	double GetMatrixtrace() const;
	void zero();
	DoubleArray DotArray(const DoubleArray& src) const;
	void Subtract(const DoubleMatrix& src, double scale = 1.);
	/*只有上半矩阵相乘: a(Tanspose) * b * scale*/
	void plusProductSymm(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1.);
	/*a(Tanspose) * b * scale*/
	void plusProductUnSymm(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1.);
	//矩阵对称化
	void symmetrized();
	//计算N矩阵
	void ComputeNMatrix(const DoubleArray& src, int dim);
	// += a * b
	void PlusProduct(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1);
	//a(tansposed) * b * scale
	void TProduct(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1);
	// this = a * b
	void beProductOf(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1);
	// this = a * b(transposed)
	void beProductTOf(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1);
	//友元重载
	friend std::ostream& operator<<(std::ostream &out, const DoubleMatrix &src);
	//将矢量扩展为矩阵
	void beMatrixForm(const DoubleArray& src);
	//在给大小之后成为单元矩阵
	void beUnitMatrix();
	/*
	 * 成为只有部分元素的矩阵，用于3d矩阵转为平面应力或平面应变矩阵
	 * a:行的索引；b :列的索引
	 */
	void beSubMatrix(const DoubleMatrix& src, const IntArray& a, const IntArray& b);
	/*
	 * 设置某矩阵为该矩阵的子矩阵, 从r行c列开始插入
	 */
	void setSubMatrix(const DoubleMatrix& src, int r, int c);
	void printfYourself() const;
	// 复制相应向量到指定的行，从相应的列开始
	void copySubVectoRow(const DoubleArray& src, int r, int c);
	// 获取一个矩阵的相应行列
	void copyMtrxOf(const DoubleMatrix& src, int rstart, int cstart, int rend, int cend);
	// 根据输入参数重设矩阵大小(这里r,c不能比原矩阵行列数目小，否则会丢失数据)
	void resizewithDatas(int r, int c);
	/*
	 * 转换矩阵： answer = r{T} * answer * r(默认是‘n’)
	 * answer = r * answer * r{T}(默认是‘t’)
	 */
	void rotatedwith(const DoubleMatrix& r, char mode = 'n');
	/*
	 * 组装矩阵
	 */
	void assemble(const DoubleMatrix& src, const IntArray& loc);
	/*
	 * 求解线性方程组 b = this * answer，使用高斯消去法
	 * 返回this（K矩阵）的奇异性，false 是K矩阵是奇异的
	 */
	bool solveForRhs(const DoubleArray& b, DoubleArray& answer, bool transpose = false);
	/**
	 * \brief 成为对角矩阵，对角元素是src中的元素
	 * \param src 对角元素
	 */
	void beDiagonal(const DoubleArray& src);
	/*
	 * 简易张量，默认是6 * 6的矩阵 
	 * follow order: 11 22 33 23 13 12
	 *  I = [1 1 1 0 0 0] X表示dyadic
	 */
	/**
	 * \brief delta_ij * delta_kl  (6 * 6) = I_ijkl---I X I
	 */
	void be_I_X_I_Matrix();
	/**
	 * \brief 0.5 * (delta_ik * delta_jl + delta_il * delta_jk) = Is_ijkl
	 */
	void beIs();
	/**
	 * \brief Id = Is - 1/3 * I X I 
	 */
	void beId();
	/**
	 * \brief a X b （向量a和b的矢量积）
	 * \param a 
	 * \param b 
	 */
	void beDyadicProductOf(const DoubleArray& a, const DoubleArray& b);
	/**
	 * \brief 叉积相加但只形成上半部矩阵
	 * \param a 向量
	 * \param scale 放大倍数
	 */
	void plusDyadicSymm(const DoubleArray& a, double scale = 1.);
	/**
	 * \brief 叉积
	 * \param a 叉乘的向量a
	 * \param b 叉乘向量b
	 * \param scale 放大倍数
	 */
	void plusDyadicUnsymm(const DoubleArray& a, double scale = 1.);

	/**
	 * \brief 成为原矩阵的对称部分 this = 0.5 * (src(i, j) + src(j, i))
	 * \param src 
	 */
	void beSymPartOf(const DoubleMatrix& src);
	/**
	 * \brief 成为原矩阵的反对称部分 this = 0.5 * (src(i, j) - src(j, i))
	 * \param src 
	 */
	void beAntiSymPartOf(const DoubleMatrix& src);
	/**
	 * \brief 计算特征向量的向量积vec X vec
	 * \param vec 特征向量
	 */
	void eigProjectionOf(const DoubleMatrix& vec);

	DoubleMatrix cutoff(DoubleMatrix A, int i, int j);

	double det(DoubleMatrix A);

	DoubleMatrix company(DoubleMatrix A);

	DoubleMatrix num_mul(DoubleMatrix A, double num);
	void matrix_inverse(DoubleMatrix& a, DoubleMatrix& b);
	DoubleMatrix inv44(DoubleMatrix t);
	DoubleMatrix inv(DoubleMatrix& a);
	void swap2(double& a, double& b);
};
	
#endif