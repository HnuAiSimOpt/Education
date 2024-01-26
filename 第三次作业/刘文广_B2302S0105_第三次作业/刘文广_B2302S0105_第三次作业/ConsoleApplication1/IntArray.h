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
	//���캯��
	IntArray(int n = 0) : values(n){}
	IntArray(const IntArray& src) : values(src.values){}
	IntArray& operator=(const IntArray& src){ values = src.values; return *this; }
	IntArray& operator=(IntArray&& src){ values = std::move(src.values); return *this; }
	IntArray(IntArray&& src) : values(std::move(src.values)) {}

	IntArray(double) = delete;
	
	//������ Ϊ��ѭ��
	std::vector<int>::iterator begin() { return this->values.begin(); }
	std::vector<int>::iterator end() { return this->values.end(); }
	std::vector<int>::const_iterator begin() const { return this->values.begin(); }
	std::vector<int>::const_iterator end() const { return this->values.end(); }
	//�б��ʼ��
	inline IntArray(std::initializer_list<int> list) : values(list) {}
	inline IntArray& operator=(std::initializer_list<int> list) { values = list; return *this; }
	
	//��������
	~IntArray(){}
	//ȡֵ
	inline int& operator()(const int& i) { return values[i - 1]; }
	inline int operator()(const int& i) const { return values[i - 1]; }
	int getSize() const { return (int)values.size(); }
	int& at(int i) { return values[i]; }
	int at(int i) const { return values[i]; }
	//��ȡֵ��ָ��͸���
	inline const int* GetPointer() const { return values.data(); }
	inline int* GetPointer() { return values.data(); }
	IntArray* GetCopy() { return new IntArray(*this); }
	//�ж��Ƿ���ڷ�����ֵ�ͷǿ�
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
	//��ֵ
	void SetValue(const int& i, const int& val) { values[i - 1] = val; }
	void push_back(const int &val) { values.push_back(val); }
	//����
	void zero(){ std::fill(this->values.begin(), this->values.end(), 0); }
	//�ı��С
	void resize(const int &i){ this->values.resize(i); }
	// ����������С
	void reserve(int futuresize) { this->values.reserve(futuresize); }
	// ���������޸Ĵ�Сn+alloc
	void resizewithvalue(int n, int alloc = 0);
	//���
	void clear() { values.clear(); }
	//����
	void sortmin2max(){ std::sort(values.begin(), values.end()); }
	void sortmax2min(){ std::sort(values.begin(), values.end(), std::greater<int>()); }
	//Ѱ�����ֵ����Сֵ
	int MaxValue() const { return *std::max_element(values.begin(), values.end()); }
	int MinValue() const { return *std::min_element(values.begin(), values.end()); }
	//ȥ����Ԫ�أ�����������Ԫ��
	void GetRidOfZero(IntArray& Input);
	// �������Ԫ�صĸ���
	int getNumberOfNonZero() const;
	//Ѱ�Ҷ�Ӧֵ������
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
	// Ѱ����Ӧֵ������ (ֵ�洢��answer��)
	void FindIndexOf(int value, IntArray& answer) const;
	//����һ�����鵽��һ��
	void FollowBy(const IntArray& src) { values.insert(values.end(), src.begin(), src.end()); }
	void FollowBy(int n, int alloc = 0);
	//�����С����
	void sort(); 
	void sort(IntArray& answer);
	//��ֵ��1��ʼ���������ֵ
	void assignfrom(int max);
	// ����ֵ����value
	void assign(int value);
	//����ĳԪ�ػ���δ��������
	bool FindValue(const int &x);
	//������ظ�����õ�����
	bool InsertSortedValueOnce(const int& x);
	//���
	void printfYouself() const;
	// ������Ԫ�����
	int addAllElements() const;
	void printfToFile(FILE* fp);
};


#endif