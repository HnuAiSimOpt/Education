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
	//����
	int nRows;
	//����
	int nColumns;
	//ֵ
	std::vector<double>values;
public:
	DoubleMatrix(int _nRows = 0, int _nColumns = 0) : nRows(_nRows), nColumns(_nColumns), values(_nRows*_nColumns)
	{
	}
	//�������캯��
	DoubleMatrix(const DoubleMatrix& mat) : nRows(mat.nRows), nColumns(mat.nColumns), values(mat.values) {}
	DoubleMatrix(DoubleMatrix&& mat) : nRows(mat.nRows), nColumns(mat.nColumns), values(std::move(mat.values)){}
	DoubleMatrix& operator=(DoubleMatrix&& mat)
	{
		nRows = mat.nRows;
		nColumns = mat.nColumns;
		values = std::move(mat.values);
		return *this;
	}
	//��ֵ����
	DoubleMatrix& operator=(const DoubleMatrix& mat)
	{
		nRows = mat.nRows;
		nColumns = mat.nColumns;
		values = mat.values;
		return *this;
	}
	~DoubleMatrix();
	//��������
	void SetValues(const int i, const int j, const double value);
	void SetnRows(const int i) { nRows = i; }
	void SetnColumns(const int j) { nColumns = j; }
	double GetValues(const int i, const int j) const { return (*this)(i, j); }
	int GetRow() const { return nRows; }
	int GetColumn() const { return nColumns; }
	void clear() { values.clear(); nRows = 0; nColumns = 0; }
	void resize(int rows, int columns) { nRows = rows; nColumns = columns; values.resize(nRows*nColumns); }
	//�б��ʼ��
	DoubleMatrix(std::initializer_list<std::initializer_list<double>>src);
	DoubleMatrix& operator=(std::initializer_list<std::initializer_list<double>>src);
	//�жϾ����Ƿ������
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
	//�ж��Ƿ�Ϊ����
	inline bool isSquare() const
	{
		return nRows == nColumns;
	}
	//�ж��Ƿ�Ϊ��
	inline bool isEmpty() const
	{
		return nRows == 0 || nColumns == 0;
	}
	//�õ�i��j��Ԫ��
	//���޸ĵ���ֵ
	inline double& operator() (const int& i, const int& j) { return values[nColumns * (i - 1) + j - 1]; }
	inline double operator() (const int& i, const int& j) const { return values[nColumns * (i - 1) + j - 1]; }
	double& at(int i, int j) { return values[nColumns * i + j]; }
	double at(int i, int j) const { return values[nColumns * i + j]; }
	//���������ʽ��ֵ
	double ComputeDeterminant() const;
	//��������
	//ע�ⲻ����ͬά��Ҳ�����
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
	/*ֻ���ϰ�������: a(Tanspose) * b * scale*/
	void plusProductSymm(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1.);
	/*a(Tanspose) * b * scale*/
	void plusProductUnSymm(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1.);
	//����Գƻ�
	void symmetrized();
	//����N����
	void ComputeNMatrix(const DoubleArray& src, int dim);
	// += a * b
	void PlusProduct(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1);
	//a(tansposed) * b * scale
	void TProduct(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1);
	// this = a * b
	void beProductOf(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1);
	// this = a * b(transposed)
	void beProductTOf(const DoubleMatrix& a, const DoubleMatrix& b, double scale = 1);
	//��Ԫ����
	friend std::ostream& operator<<(std::ostream &out, const DoubleMatrix &src);
	//��ʸ����չΪ����
	void beMatrixForm(const DoubleArray& src);
	//�ڸ���С֮���Ϊ��Ԫ����
	void beUnitMatrix();
	/*
	 * ��Ϊֻ�в���Ԫ�صľ�������3d����תΪƽ��Ӧ����ƽ��Ӧ�����
	 * a:�е�������b :�е�����
	 */
	void beSubMatrix(const DoubleMatrix& src, const IntArray& a, const IntArray& b);
	/*
	 * ����ĳ����Ϊ�þ�����Ӿ���, ��r��c�п�ʼ����
	 */
	void setSubMatrix(const DoubleMatrix& src, int r, int c);
	void printfYourself() const;
	// ������Ӧ������ָ�����У�����Ӧ���п�ʼ
	void copySubVectoRow(const DoubleArray& src, int r, int c);
	// ��ȡһ���������Ӧ����
	void copyMtrxOf(const DoubleMatrix& src, int rstart, int cstart, int rend, int cend);
	// �������������������С(����r,c���ܱ�ԭ����������ĿС������ᶪʧ����)
	void resizewithDatas(int r, int c);
	/*
	 * ת������ answer = r{T} * answer * r(Ĭ���ǡ�n��)
	 * answer = r * answer * r{T}(Ĭ���ǡ�t��)
	 */
	void rotatedwith(const DoubleMatrix& r, char mode = 'n');
	/*
	 * ��װ����
	 */
	void assemble(const DoubleMatrix& src, const IntArray& loc);
	/*
	 * ������Է����� b = this * answer��ʹ�ø�˹��ȥ��
	 * ����this��K���󣩵������ԣ�false ��K�����������
	 */
	bool solveForRhs(const DoubleArray& b, DoubleArray& answer, bool transpose = false);
	/**
	 * \brief ��Ϊ�ԽǾ��󣬶Խ�Ԫ����src�е�Ԫ��
	 * \param src �Խ�Ԫ��
	 */
	void beDiagonal(const DoubleArray& src);
	/*
	 * ����������Ĭ����6 * 6�ľ��� 
	 * follow order: 11 22 33 23 13 12
	 *  I = [1 1 1 0 0 0] X��ʾdyadic
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
	 * \brief a X b ������a��b��ʸ������
	 * \param a 
	 * \param b 
	 */
	void beDyadicProductOf(const DoubleArray& a, const DoubleArray& b);
	/**
	 * \brief �����ӵ�ֻ�γ��ϰ벿����
	 * \param a ����
	 * \param scale �Ŵ���
	 */
	void plusDyadicSymm(const DoubleArray& a, double scale = 1.);
	/**
	 * \brief ���
	 * \param a ��˵�����a
	 * \param b �������b
	 * \param scale �Ŵ���
	 */
	void plusDyadicUnsymm(const DoubleArray& a, double scale = 1.);

	/**
	 * \brief ��Ϊԭ����ĶԳƲ��� this = 0.5 * (src(i, j) + src(j, i))
	 * \param src 
	 */
	void beSymPartOf(const DoubleMatrix& src);
	/**
	 * \brief ��Ϊԭ����ķ��ԳƲ��� this = 0.5 * (src(i, j) - src(j, i))
	 * \param src 
	 */
	void beAntiSymPartOf(const DoubleMatrix& src);
	/**
	 * \brief ��������������������vec X vec
	 * \param vec ��������
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