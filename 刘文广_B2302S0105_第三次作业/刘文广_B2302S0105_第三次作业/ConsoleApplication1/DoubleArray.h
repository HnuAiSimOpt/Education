#ifndef DOUBLEARRAY_H_
#define DOUBLEARRAY_H_
#include <vector>
#include <cmath>
#include <numeric>
#include "InfoOut.h"
class IntArray;
class DoubleMatrix;
/*
 * 排列方式是按行排列
 */
class DoubleArray
{
protected:
	std::vector<double> values;
public:
	//构造函数
	DoubleArray(int n = 0) : values(n){}
	DoubleArray(const DoubleArray& src) : values(src.values){}
	DoubleArray& operator=(const DoubleArray& src) { values = src.values; return *this; }
	DoubleArray(DoubleArray&& src) : values(std::move(src.values)){}
	DoubleArray& operator=(DoubleArray&& src) { values = std::move(src.values); return *this; }
	//禁用的构造函数
	DoubleArray(double) = delete;
	//析构函数
	~DoubleArray(){}
	//迭代器 为了循环
	std::vector<double>::iterator begin() { return this->values.begin(); }
	std::vector<double>::iterator end() { return this->values.end(); }
	std::vector<double>::const_iterator begin() const { return this->values.begin(); }
	std::vector<double>::const_iterator end() const { return this->values.end(); }
	//获取指针和复制
	//获取值的指针和复制
	inline const double* GetPointer() const { return values.data(); }
	inline double* GetPointer() { return values.data(); }
	DoubleArray* GetCopy() { return new DoubleArray(*this); }
	//取值
	inline double& operator()(const int& i) { return values[i - 1]; }
	inline double operator() (const int& i) const { return values[i - 1]; }
	int GetSize() const { return (int)values.size();}
	double& at(int i) { return values[i]; }
	double at(int i) const { return values[i]; }
	//列表初始化
	inline DoubleArray(std::initializer_list<double> list) : values(list) {}
	inline DoubleArray& operator=(std::initializer_list<double> list) { values = list; return *this; }
	//判断是否存在非有限值和非空
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
	//赋值
	void SetValue(const int& i, const double& val) { values[i - 1] = val; }
	void push_back(const double &val) { values.push_back(val); }
	void reserve(int s);
	//重设大小
	void resize(const int s);
	// 重新设置大小，如果大小不匹配就会相应修改
	void resizewithvalue(int n, int alloc = 0);
	//置零
	void zero();
	//清空
	void clear() { values.clear(); }
	//向量运算
	//加法
	// 每个元素都加上相同的值
	void Add(const double &val);
	void Add(const DoubleArray &src);
	void Add(const int& i, const double &val);
	void Add(const double& factor, const DoubleArray& b);
	double addAllElements() const;
	DoubleArray operator+ (const DoubleArray&);
	DoubleArray operator+ (const double &);
	// this += a(T) * b * s 
	void plus_tproduct(const DoubleMatrix& a, const DoubleArray& b, double s = 1.);
	// 赋值全部为相同的值
	void assign(double value);
	//减法
	void minus(const double &val);
	void minus(const DoubleArray &src);
	void minus(const int& i, const double& val);
	DoubleArray operator- (const double&);
	DoubleArray operator- (const DoubleArray&) const;
	void Subtract(const DoubleArray& a, const DoubleArray& b);
	void Subtract(const DoubleArray& a);
	void beDifferenceOf(const DoubleArray& a, const DoubleArray& b);
	// 幂运算
	void power(const double& x);
	//点乘
	double Dot(const DoubleArray&) const;
	// 在两个向量中选择最小或最大的
	void beMinOf(const DoubleArray& a, const DoubleArray& b);
	void beMaxOf(const DoubleArray& a, const DoubleArray& b);
	//两向量想乘得出一个矩阵
	DoubleMatrix operator*(const DoubleArray&);
	//两向量叉积 this = a X v(注意这里哪个在前面哪个在后面)
	void VectorProduct(const DoubleArray& a, const DoubleArray &v);
	// 复制一个向量中的元素到末尾
	void copySubVec(const DoubleArray& src, int si);
	//Norm
	double SquareNorm();
	double SquareNorm() const;
	// 向量平方范数开方
	double Norm();
	//单位化
	double Normlize();
	//扩大一定的倍数
	void Scale(double s);
	// this = src * s
	void beScaledOf(const DoubleArray& src, double s);
	//矩阵 * 向量 = this向量(a需要转置)
	void MtrxTDotVec(const DoubleMatrix& a, const DoubleArray& b);
	void MtrxDotVect(const DoubleMatrix& a, const DoubleArray& b);
	//组装
	void assemble(const DoubleArray& a, const IntArray& loc);
	// 组装向量的平方
	void assembleSquared(const DoubleArray& a, const IntArray& loc);
	//插入一个数组到另一个
	void FollowBy(const DoubleArray& src) { values.insert(values.end(), src.begin(), src.end()); }
	//两向量之间的距离
	double Length(const DoubleArray& src) const;
	double Length(const DoubleArray* src) const;
	static double LengthCoords(double x1, double x2, double y1, double y2, double z1 = 0, double z2 = 0);
	// 将对称应力矩阵转换为向量
	void beSymVectorOfStress(const DoubleMatrix& src);
	//将对称应变矩阵转换为向量
	void beSymVectorOfStrain(const DoubleMatrix& src);
	//将矩阵转换为向量
	void beVectorForm(const DoubleMatrix& src);
	/* 
	 * 对向量施加旋转，默认是n 
	 * this = r{T} * this(mode = 'n'); this = r * this(mode = 't')
	 * 局部 = r * 全局坐标 --- mode = 't'
	 * 全局坐标 = r{T} * 局部坐标  --- mode = 'n'
	 */
	void rotatewith(const DoubleMatrix& r, char mode = 'n');
	// 从矩阵获取某一行向量(r行cstart列开始到cend列结束)
	void copyRowVectorOf(const DoubleMatrix& src, int r, int cstart, int cend);
	// 从矩阵获取某一列向量(c列rstart行开始到rend行结束)
	void copyColumnVectorOf(const DoubleMatrix& src, int c, int rstart, int rend);
	//输出
	void printYourself() const;
	// 把所有元素置负
	void negated();
	// 特征值专用，将omega^2转换为频率
	void getFrequentcy(DoubleArray& fre);
	// 判断所有元素是不是全为0
	bool containsOnlyZeros() const;
	// 输出到文件
	void printfToFile(FILE* fp);


	// double = a(T) * M * a
	double calenergy(const DoubleMatrix& src);
};


#endif