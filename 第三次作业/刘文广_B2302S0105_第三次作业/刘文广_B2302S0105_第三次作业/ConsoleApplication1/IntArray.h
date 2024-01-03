#ifndef INTARRAY_H_
#define INTARRAY_H_
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
using namespace std;

class IntArray
{
private:
	std::vector<int> values;
public:
	//构造函数
	IntArray(int n = 0) : values(n){}
	IntArray(const IntArray& src) : values(src.values){}
	IntArray& operator=(const IntArray& src){ values = src.values; return *this; }
	IntArray& operator=(IntArray&& src){ values = std::move(src.values); return *this; }
	IntArray(IntArray&& src) : values(std::move(src.values)) {}

	IntArray(double) = delete;
	
	//迭代器 为了循环
	std::vector<int>::iterator begin() { return this->values.begin(); }
	std::vector<int>::iterator end() { return this->values.end(); }
	std::vector<int>::const_iterator begin() const { return this->values.begin(); }
	std::vector<int>::const_iterator end() const { return this->values.end(); }
	//列表初始化
	inline IntArray(std::initializer_list<int> list) : values(list) {}
	inline IntArray& operator=(std::initializer_list<int> list) { values = list; return *this; }
	
	//析构函数
	~IntArray(){}
	//取值
	inline int& operator()(const int& i) { return values[i - 1]; }
	inline int operator()(const int& i) const { return values[i - 1]; }
	int getSize() const { return (int)values.size(); }
	int& at(int i) { return values[i]; }
	int at(int i) const { return values[i]; }
	//获取值的指针和复制
	inline const int* GetPointer() const { return values.data(); }
	inline int* GetPointer() { return values.data(); }
	IntArray* GetCopy() { return new IntArray(*this); }
	//判断是否存在非有限值和非空
	inline bool isFinite() const
	{
		for (int val : values)
		{
			if (!std::isfinite((double)val))
			{
				return false;
			}
		}
		return true;
	}
	inline bool isEmpty() const { return this->values.empty(); }
	//赋值
	void SetValue(const int& i, const int& val) { values[i - 1] = val; }
	void push_back(const int &val) { values.push_back(val); }
	//置零
	void zero(){ std::fill(this->values.begin(), this->values.end(), 0); }
	//改变大小
	void resize(const int &i){ this->values.resize(i); }
	// 设置容量大小
	void reserve(int futuresize) { this->values.reserve(futuresize); }
	// 根据输入修改大小n+alloc
	void resizewithvalue(int n, int alloc = 0);
	//清除
	void clear() { values.clear(); }
	//排序
	void sortmin2max(){ std::sort(values.begin(), values.end()); }
	void sortmax2min(){ std::sort(values.begin(), values.end(), std::greater<int>()); }
	//寻找最大值和最小值
	int MaxValue() const { return *std::max_element(values.begin(), values.end()); }
	int MinValue() const { return *std::min_element(values.begin(), values.end()); }
	//去掉零元素，并保留非零元素
	void GetRidOfZero(IntArray& Input);
	// 计算非零元素的个数
	int getNumberOfNonZero() const;
	//寻找对应值的索引
	int findFirstIndexOf(int value) const
	{
		auto it = std::find(values.begin(), values.end(), value);
		if (it == values.end())
		{
			return 0;
		}
		else
		{
			return (int)(it - values.begin() + 1);
		}
	}
	bool contains(int val) const { return findFirstIndexOf(val) > 0; }
	// 寻找相应值的索引 (值存储在answer中)
	void FindIndexOf(int value, IntArray& answer) const;
	//插入一个数组到另一个
	void FollowBy(const IntArray& src) { values.insert(values.end(), src.begin(), src.end()); }
	void FollowBy(int n, int alloc = 0);
	//排序从小到大
	void sort(); 
	void sort(IntArray& answer);
	//赋值从1开始到输入最大值
	void assignfrom(int max);
	// 所有值赋予value
	void assign(int value);
	//查找某元素基于未排序数组
	bool FindValue(const int &x);
	//插入非重复排序好的数组
	bool InsertSortedValueOnce(const int& x);
	//输出
	void printfYouself() const;
	// 把所有元素相加
	int addAllElements() const;
	void printfToFile(FILE* fp);
};


#endif