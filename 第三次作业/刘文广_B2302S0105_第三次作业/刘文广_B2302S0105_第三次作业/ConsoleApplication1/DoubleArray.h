#ifndef DOUBLEARRAY_H_
#define DOUBLEARRAY_H_
#include <vector>
#include <cmath>
#include <numeric>
#include "InfoOut.h"
class IntArray;
class DoubleMatrix;
/*
 * ���з�ʽ�ǰ�������
 */
class DoubleArray
{
protected:
	std::vector<double> values;
public:
	//���캯��
	DoubleArray(int n = 0) : values(n){}
	DoubleArray(const DoubleArray& src) : values(src.values){}
	DoubleArray& operator=(const DoubleArray& src) { values = src.values; return *this; }
	DoubleArray(DoubleArray&& src) : values(std::move(src.values)){}
	DoubleArray& operator=(DoubleArray&& src) { values = std::move(src.values); return *this; }
	//���õĹ��캯��
	DoubleArray(double) = delete;
	//��������
	~DoubleArray(){}
	//������ Ϊ��ѭ��
	std::vector<double>::iterator begin() { return this->values.begin(); }
	std::vector<double>::iterator end() { return this->values.end(); }
	std::vector<double>::const_iterator begin() const { return this->values.begin(); }
	std::vector<double>::const_iterator end() const { return this->values.end(); }
	//��ȡָ��͸���
	//��ȡֵ��ָ��͸���
	inline const double* GetPointer() const { return values.data(); }
	inline double* GetPointer() { return values.data(); }
	DoubleArray* GetCopy() { return new DoubleArray(*this); }
	//ȡֵ
	inline double& operator()(const int& i) { return values[i - 1]; }
	inline double operator() (const int& i) const { return values[i - 1]; }
	int GetSize() const { return (int)values.size();}
	double& at(int i) { return values[i]; }
	double at(int i) const { return values[i]; }
	//�б��ʼ��
	inline DoubleArray(std::initializer_list<double> list) : values(list) {}
	inline DoubleArray& operator=(std::initializer_list<double> list) { values = list; return *this; }
	//�ж��Ƿ���ڷ�����ֵ�ͷǿ�
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
	inline bool isEmpty() const { return this->values.empty(); }
	//��ֵ
	void SetValue(const int& i, const double& val) { values[i - 1] = val; }
	void push_back(const double &val) { values.push_back(val); }
	void reserve(int s);
	//�����С
	void resize(const int s);
	// �������ô�С�������С��ƥ��ͻ���Ӧ�޸�
	void resizewithvalue(int n, int alloc = 0);
	//����
	void zero();
	//���
	void clear() { values.clear(); }
	//��������
	//�ӷ�
	// ÿ��Ԫ�ض�������ͬ��ֵ
	void Add(const double &val);
	void Add(const DoubleArray &src);
	void Add(const int& i, const double &val);
	void Add(const double& factor, const DoubleArray& b);
	double addAllElements() const;
	DoubleArray operator+ (const DoubleArray&);
	DoubleArray operator+ (const double &);
	// this += a(T) * b * s 
	void plus_tproduct(const DoubleMatrix& a, const DoubleArray& b, double s = 1.);
	// ��ֵȫ��Ϊ��ͬ��ֵ
	void assign(double value);
	//����
	void minus(const double &val);
	void minus(const DoubleArray &src);
	void minus(const int& i, const double& val);
	DoubleArray operator- (const double&);
	DoubleArray operator- (const DoubleArray&) const;
	void Subtract(const DoubleArray& a, const DoubleArray& b);
	void Subtract(const DoubleArray& a);
	void beDifferenceOf(const DoubleArray& a, const DoubleArray& b);
	// ������
	void power(const double& x);
	//���
	double Dot(const DoubleArray&) const;
	// ������������ѡ����С������
	void beMinOf(const DoubleArray& a, const DoubleArray& b);
	void beMaxOf(const DoubleArray& a, const DoubleArray& b);
	//��������˵ó�һ������
	DoubleMatrix operator*(const DoubleArray&);
	//��������� this = a X v(ע�������ĸ���ǰ���ĸ��ں���)
	void VectorProduct(const DoubleArray& a, const DoubleArray &v);
	// ����һ�������е�Ԫ�ص�ĩβ
	void copySubVec(const DoubleArray& src, int si);
	//Norm
	double SquareNorm();
	double SquareNorm() const;
	// ����ƽ����������
	double Norm();
	//��λ��
	double Normlize();
	//����һ���ı���
	void Scale(double s);
	// this = src * s
	void beScaledOf(const DoubleArray& src, double s);
	//���� * ���� = this����(a��Ҫת��)
	void MtrxTDotVec(const DoubleMatrix& a, const DoubleArray& b);
	void MtrxDotVect(const DoubleMatrix& a, const DoubleArray& b);
	//��װ
	void assemble(const DoubleArray& a, const IntArray& loc);
	// ��װ������ƽ��
	void assembleSquared(const DoubleArray& a, const IntArray& loc);
	//����һ�����鵽��һ��
	void FollowBy(const DoubleArray& src) { values.insert(values.end(), src.begin(), src.end()); }
	//������֮��ľ���
	double Length(const DoubleArray& src) const;
	double Length(const DoubleArray* src) const;
	static double LengthCoords(double x1, double x2, double y1, double y2, double z1 = 0, double z2 = 0);
	// ���Գ�Ӧ������ת��Ϊ����
	void beSymVectorOfStress(const DoubleMatrix& src);
	//���Գ�Ӧ�����ת��Ϊ����
	void beSymVectorOfStrain(const DoubleMatrix& src);
	//������ת��Ϊ����
	void beVectorForm(const DoubleMatrix& src);
	/* 
	 * ������ʩ����ת��Ĭ����n 
	 * this = r{T} * this(mode = 'n'); this = r * this(mode = 't')
	 * �ֲ� = r * ȫ������ --- mode = 't'
	 * ȫ������ = r{T} * �ֲ�����  --- mode = 'n'
	 */
	void rotatewith(const DoubleMatrix& r, char mode = 'n');
	// �Ӿ����ȡĳһ������(r��cstart�п�ʼ��cend�н���)
	void copyRowVectorOf(const DoubleMatrix& src, int r, int cstart, int cend);
	// �Ӿ����ȡĳһ������(c��rstart�п�ʼ��rend�н���)
	void copyColumnVectorOf(const DoubleMatrix& src, int c, int rstart, int rend);
	//���
	void printYourself() const;
	// ������Ԫ���ø�
	void negated();
	// ����ֵר�ã���omega^2ת��ΪƵ��
	void getFrequentcy(DoubleArray& fre);
	// �ж�����Ԫ���ǲ���ȫΪ0
	bool containsOnlyZeros() const;
	// ������ļ�
	void printfToFile(FILE* fp);


	// double = a(T) * M * a
	double calenergy(const DoubleMatrix& src);
};


#endif