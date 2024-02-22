#include "IntArray.h"

void IntArray::FollowBy(int n, int alloc)
{
	if (alloc && (int)values.capacity() < this->getSize() + 1)
	{
		values.reserve(values.capacity() + alloc + 1);
	}
	values.push_back(n);
}

void IntArray::sort()
{
	std::sort(this->begin(), this->end());
}
void IntArray::sort(IntArray& answer)
{
	answer.clear();
	answer = *this;
	answer.sort();
}
bool IntArray::FindValue(const int& x)
{
	IntArray tmp;
	this->sort(tmp);
	return std::binary_search(tmp.begin(), tmp.end(), x) == 1 ? true : false;
}
bool IntArray::InsertSortedValueOnce(const int& x)
{
	auto low = std::lower_bound(values.begin(), values.end(), x);
	if (low == values.end() || *low != x)
	{
		values.insert(low, x);
		return true;
	}
	return false;
}

void IntArray::assignfrom(int max)
{
	this->resize(max);
	for (int i = 1; i <= max; i++)
	{
		(*this)(i) = i;
	}
}

void IntArray::assign(int value)
{
	for (int i = 1; i <= values.size(); i++)
	{
		(*this)(i) = value;
	}
}

void IntArray::printfYouself() const
{
	printf("size: %d\n", this->getSize());
	for (int x : *this)
	{
		printf("%5d", x);
	}
	printf("\n");
}

void IntArray::resizewithvalue(int n, int alloc)
{
	if (alloc > 0 && (int)this->values.capacity() < n)
	{
		this->values.reserve(n + alloc);
	}
	this->values.resize(n);
}

int IntArray::addAllElements() const
{
	int sum = 0;
	for (int i = 0; i < values.size(); i++)
	{
		sum += values[i];
	}
	return sum;
}

void IntArray::printfToFile(FILE* fp)
{
	int m = values.size() / 16;
	int rest = values.size() - m * 16;
	for (int i = 1; i <= m; i++)
	{
		for (int j = 1; j <= 16; j++)
		{
			fprintf(fp, "%16d", (*this)((i - 1) * 16 + j));
		}
		fprintf(fp, "\n");
	}
	for (int i = 1; i <= rest; i++)
	{
		fprintf(fp, "%16d", (*this)(16 * m + i));
	}
	fprintf(fp, "\n");
}

void IntArray::GetRidOfZero(IntArray& Input)
{
	int size = 0;
	for (const int& x : Input.values)
	{
		if (x)
		{
			size++;
		}
	}
	int pos = 1;
	this->values.resize(size);
	for (int i = 1; i <= Input.getSize(); i++)
	{
		if (Input(i))
		{
			(*this)(pos++) = i;
		}
	}
}

int IntArray::getNumberOfNonZero() const
{
	int size = 0;
	for (const int& x : this->values)
	{
		if (x)
		{
			size++;
		}
	}
	return size;
}


void IntArray::FindIndexOf(int value, IntArray& answer) const
{
	for (int i = 1; i <= values.size(); i++)
	{
		if (value == values[i - 1])
		{
			answer.push_back(i);
		}
	}
}
