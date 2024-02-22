#include "DoubleMatrix.h"
#include "DoubleArray.h"
#include "IntArray.h"
#include "mathfem.h"
#include <iostream>
DoubleMatrix::~DoubleMatrix(){}
double DoubleMatrix::ComputeDeterminant() const
{
	if (!this->isSquare())
	{
		ERROR("Not square matrix");
	}
	if (nRows == 1)
	{
		return values[0];
	}
	else if (nRows == 2)
	{
		return (values[0] * values[3] - values[1] * values[2]);
	}
	else if (nRows == 3)
	{
		return (values[0] * values[4] * values[8] + values[3] * values[7] * values[2]
			  + values[1] * values[5] * values[6] - values[2] * values[4] * values[6]
			  - values[1] * values[3] * values[8] - values[0] * values[5] * values[7]);
	}
	else
	{
		ERROR("cannot compute the determinant of a matrix larger than 3x3");
	}
	//return 0.0;
}
DoubleMatrix DoubleMatrix::operator+(const DoubleMatrix& src) const
{
	if (nRows != src.nRows||nColumns != src.nColumns)
	{
		ERROR("cannot add the matrixs with different row or columns");
	}
	DoubleMatrix tmp(src.nRows, src.nColumns);
	for (int i = 1; i <= src.nRows; i++)
	{
		for (int j = 1; j <= src.nColumns; j++)
		{
			tmp(i, j) = (*this)(i, j) + src(i, j);
		}
	}
	return tmp;

}
DoubleMatrix DoubleMatrix::operator-(const DoubleMatrix &src) const
{
	if (nRows != src.nRows || nColumns != src.nColumns)
	{
		ERROR("cannot minus the matrixs with different row or columns");
	}
	DoubleMatrix tmp(src.nRows, src.nColumns);
	for (int i = 1; i <= src.nRows; i++)
	{
		for (int j = 1; j <= src.nColumns; j++)
		{
			tmp(i, j) = (*this)(i, j) - src(i, j);
		}
	}
	return tmp;
}

DoubleMatrix DoubleMatrix::operator*(const DoubleMatrix &src) const
{
	if (nColumns != src.nRows)
	{
		ERROR("cannot mutiply the matrixs with different row or columns");
	}
	DoubleMatrix tmp(nRows, src.nColumns);
	tmp.zero();
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= src.nColumns; j++)
		{
			for (int k = 1; k <= nColumns; k++)
			{
				tmp(i, j) += (*this)(i, k) * src(k, j);
			}
			
		}
	}
	return tmp;
}
DoubleMatrix DoubleMatrix::operator*(const double& x) const
{
	DoubleMatrix tmp(nRows, nColumns);
	tmp.zero();
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			tmp(i, j) = (*this)(i, j)*x;
		}
	}
	return tmp;
}
void DoubleMatrix::zero()
{
	std::fill(this->values.begin(), this->values.end(), 0.0);
}
DoubleMatrix DoubleMatrix::Inverse() const
{
	DoubleMatrix tmp(nRows, nColumns);
	if (!isSquare())
	{
		ERROR("cannot inverse");
	}
	if (nRows == 2)
	{
		double det = this->ComputeDeterminant();
		if (det == 0.0)
		{
			ERROR("det==0");
		}
		tmp(1, 1) = values[3] / det;
		tmp(1, 2) = -values[1] / det;
		tmp(2, 1) = -values[2] / det;
		tmp(2, 2) = values[0] / det;
		return tmp;
	}
	else if (nRows == 3)
	{
		
		double det = this->ComputeDeterminant();
		if (det == 0.0)
		{
			ERROR("det==0");
		}
		tmp(1, 1) = (values[4] * values[8] - values[5] * values[7]) / det;
		tmp(1, 2) = -(values[1] * values[8] - values[7] * values[2]) / det;
		tmp(1, 3) = (values[1] * values[5] - values[4] * values[2]) / det;
		tmp(2, 1) = -(values[3] * values[8] - values[5] * values[6]) / det;
		tmp(2, 2) = (values[0] * values[8] - values[6] * values[2]) / det;
		tmp(2, 3) = -(values[0] * values[5] - values[3] * values[2]) / det;
		tmp(3, 1) = (values[3] * values[7] - values[4] * values[6]) / det;
		tmp(3, 2) = -(values[0] * values[7] - values[1] * values[6]) / det;
		tmp(3, 3) = (values[0] * values[4] - values[1] * values[3]) / det;
		return tmp;
	}
	return tmp;
}

void DoubleMatrix::beInverseOf(const DoubleMatrix& src)
{
	double det;
	if (!src.isSquare())
	{
		ERROR("cannot inverse matrix with different columns and rows");
	}
	this->resize(src.nRows, src.nColumns);
	if (nRows == 1)
	{
		(*this)(1, 1) = 1. / src(1, 1);
	}
	else if (nRows == 2)
	{
		det = src(1, 1) * src(2, 2) - src(1, 2) * src(2, 1);
		if (det == 0.0)
		{
			ERROR("det==0");
		}
		(*this)(1, 1) = src(2, 2) / det;
		(*this)(1, 2) = -src(1, 2) / det;
		(*this)(2, 1) = -src(2, 1)/ det;
		(*this)(2, 2) = src(1, 1) / det;
	}
	else if (nRows == 3)
	{
		det = src(1, 1) * src(2, 2) * src(3, 3) + src(1, 2) * src(2, 3) * src(3, 1) +
			src(1, 3) * src(2, 1) * src(3, 2) - src(1, 3) * src(2, 2) * src(3, 1) -
			src(2, 3) * src(3, 2) * src(1, 1) - src(3, 3) * src(1, 2) * src(2, 1);

		(*this)(1, 1) = (src(2, 2) * src(3, 3) - src(2, 3) * src(3, 2)) / det;
		(*this)(2, 1) = (src(2, 3) * src(3, 1) - src(2, 1) * src(3, 3)) / det;
		(*this)(3, 1) = (src(2, 1) * src(3, 2) - src(2, 2) * src(3, 1)) / det;
		(*this)(1, 2) = (src(1, 3) * src(3, 2) - src(1, 2) * src(3, 3)) / det;
		(*this)(2, 2) = (src(1, 1) * src(3, 3) - src(1, 3) * src(3, 1)) / det;
		(*this)(3, 2) = (src(1, 2) * src(3, 1) - src(1, 1) * src(3, 2)) / det;
		(*this)(1, 3) = (src(1, 2) * src(2, 3) - src(1, 3) * src(2, 2)) / det;
		(*this)(2, 3) = (src(1, 3) * src(2, 1) - src(1, 1) * src(2, 3)) / det;
		(*this)(3, 3) = (src(1, 1) * src(2, 2) - src(1, 2) * src(2, 1)) / det;
	}
	else
	{
		// size >3 ... gaussian elimination - slow but safe
		//
		double piv, linkomb;
		DoubleMatrix tmp = src;
		// initialize answer to be unity matrix;
		this->zero();
		for (int i = 1; i <= nRows; i++) 
		{
			(*this)(i, i) = 1.0;
		}

		// lower triangle elimination by columns
		for (int i = 1; i < nRows; i++) 
		{

			piv = tmp(i, i);
			if (fabs(piv) < 1.e-24) 
			{
				ERROR("pivot (%d,%d) to close to small (< 1.e-24)", i, i);
			}

			for (int j = i + 1; j <= nRows; j++) 
			{
				linkomb = tmp(j, i) / tmp(i, i);
				for (int k = i; k <= nRows; k++) 
				{
					tmp(j, k) -= tmp(i, k) * linkomb;
				}

				for (int k = 1; k <= nRows; k++) {
					(*this)(j, k) -= (*this)(i, k) * linkomb;
				}
			}
		}

		// upper triangle elimination by columns
		for (int i = nRows; i > 1; i--) 
		{
			piv = tmp(i, i);
			for (int j = i - 1; j > 0; j--) 
			{
				linkomb = tmp(j, i) / piv;
				for (int k = i; k > 0; k--) 
				{
					tmp(j, k) -= tmp(i, k) * linkomb;
				}

				for (int k = nRows; k > 0; k--) 
				{
					// tmp -> at(j,k)-= tmp  ->at(i,k)*linkomb;
					(*this)(j, k) -= (*this)(i, k) * linkomb;
				}
			}
		}

		// diagonal scaling
		for (int i = 1; i <= nRows; i++) 
		{
			for (int j = 1; j <= nRows; j++) {
				(*this)(i, j) /= tmp(i, i);
			}
		}
	}

}

void DoubleMatrix::beInverseAndTransposeOf(const DoubleMatrix& src)
{
	double det;
	if (!src.isSquare())
	{
		ERROR("cannot inverse matrix with different columns and rows");
	}
	this->resize(src.nRows, src.nColumns);
	if (nRows == 1)
	{
		(*this)(1, 1) = 1. / src(1, 1);
	}
	else if (nRows == 2)
	{
		det = src(1, 1) * src(2, 2) - src(1, 2) * src(2, 1);
		if (det == 0.0)
		{
			ERROR("det==0");
		}
		(*this)(1, 1) = src(2, 2) / det;
		(*this)(1, 2) = -src(2, 1) / det;
		(*this)(2, 1) = -src(1, 2) / det;
		(*this)(2, 2) = src(1, 1) / det;
	}
	else if (nRows == 3)
	{
		det = src(1, 1) * src(2, 2) * src(3, 3) + src(1, 2) * src(2, 3) * src(3, 1) +
			src(1, 3) * src(2, 1) * src(3, 2) - src(1, 3) * src(2, 2) * src(3, 1) -
			src(2, 3) * src(3, 2) * src(1, 1) - src(3, 3) * src(1, 2) * src(2, 1);

		(*this)(1, 1) = (src(2, 2) * src(3, 3) - src(2, 3) * src(3, 2)) / det;
		(*this)(1, 2) = (src(2, 3) * src(3, 1) - src(2, 1) * src(3, 3)) / det;
		(*this)(1, 3) = (src(2, 1) * src(3, 2) - src(2, 2) * src(3, 1)) / det;
		(*this)(2, 1) = (src(1, 3) * src(3, 2) - src(1, 2) * src(3, 3)) / det;
		(*this)(2, 2) = (src(1, 1) * src(3, 3) - src(1, 3) * src(3, 1)) / det;
		(*this)(2, 3) = (src(1, 2) * src(3, 1) - src(1, 1) * src(3, 2)) / det;
		(*this)(3, 1) = (src(1, 2) * src(2, 3) - src(1, 3) * src(2, 2)) / det;
		(*this)(3, 2) = (src(1, 3) * src(2, 1) - src(1, 1) * src(2, 3)) / det;
		(*this)(3, 3) = (src(1, 1) * src(2, 2) - src(1, 2) * src(2, 1)) / det;
	}
	else
	{
		ERROR("exceed the size 3!");
	}
}

DoubleMatrix DoubleMatrix::Tanspose() const
{
	DoubleMatrix tmp(nColumns, nRows);
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			tmp(j, i) = (*this)(i, j);
		}
	}
	return tmp;
}
double DoubleMatrix::GetMatrixtrace() const
{
	if (!this->isSquare())
	{
		ERROR("Not square matrix");
	}
	double answer = 0.0;
	for (int i = 1; i <= nRows; i++)
	{
		answer += (*this)(i,i);
	}
	return answer;
}
std::ostream& operator<<(std::ostream &out, const DoubleMatrix &src)
{
	for (int i = 1; i <= src.nRows; i++)
	{
		for (int j = 1; j <= src.nColumns; j++)
		{
			//cout << " " << src(i,j);
		}
		//cout << "\n";
	}
	return out;

}
void DoubleMatrix::SetValues(const int i, const int j, const double value)
{
	(*this)(i, j) = value;
}
void DoubleMatrix::Add(const DoubleMatrix& src, double s)
{
	if (this->isEmpty())
	{
		this->resize(src.nRows, src.nColumns);
		this->zero();
	}
	if (nRows < src.nRows || nColumns < src.nColumns)
	{
		ERROR("cannot add the matrixs with large row or columns");
	}
	for (int i = 1; i <= src.nRows; i++)
	{
		for (int j = 1; j <= src.nColumns; j++)
		{
			(*this)(i, j) += src(i, j) * s;
		}
	}
}
void DoubleMatrix::minus(const DoubleMatrix& src)
{
	if (nRows != src.nRows || nColumns != src.nColumns)
	{
		ERROR("cannot minus the matrixs with different row or columns");
	}
	for (int i = 1; i <=nRows; i++)
	{
		for (int j = 1; j <=nColumns; j++)
		{
			(*this)(i, j) -= src(i, j);
		}
	}
}
void DoubleMatrix::Scale(const double& cof)
{
	for (auto& val : values)
	{
		val *= cof;
	}
}

DoubleArray DoubleMatrix::DotArray(const DoubleArray& src) const
{
	if (this->nColumns != src.GetSize())
	{
		ERROR("Cannot Dot!");
	}
	DoubleArray answer(this->nRows);
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			answer(i) += (*this)(i, j)*src(j);
		}
	}
	return answer;
}

void DoubleMatrix::Subtract(const DoubleMatrix& src, double scale)
{
	if (nRows != src.nRows && nColumns != src.nColumns)
	{
		ERROR("cannot subtract for mismatch size!");
	}
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			(*this)(i, j) -= src(i, j) * scale;
		}
	}
}

void DoubleMatrix::ComputeNMatrix(const DoubleArray& src, int dim)
{
	this->resize(dim, dim * src.GetSize());
	for (int i = 1; i <= src.GetSize(); i++)
	{
		for (int j = 1; j <= dim; j++)
		{
			(*this)(j, (i - 1) * dim + j) = src(i);
		}
	}
}

void DoubleMatrix::TProduct(const DoubleMatrix& a, const DoubleMatrix& b, double scale)
{
	this->resize(a.GetColumn(), b.GetColumn());
	for (int i = 1; i <= a.GetColumn(); i++)
	{
		for (int j = 1; j <= b.GetColumn(); j++)
		{
			double tmp = 0.;
			for (int k = 1; k <= a.GetRow(); k++)
			{
				tmp += a(k, i) * b(k, j);
			}
			(*this)(i, j) = tmp;
		}
	}
	this->Scale(scale);
}

void DoubleMatrix::beProductOf(const DoubleMatrix& a, const DoubleMatrix& b, double scale)
{
	if (this->isEmpty())
	{
		this->nRows = a.GetRow();
		this->nColumns = b.GetColumn();
		this->values.assign(a.GetRow() * b.GetColumn(), 0);
	}
	this->clear();
	this->resize(a.GetRow(), b.GetColumn());
	this->zero();
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			double sum = 0.;
			for (int k = 1; k <= a.GetColumn(); k++)
			{
				sum += a(i, k) * b(k, j);
			}
			(*this)(i, j) = sum * scale;
		}
	}
}

void DoubleMatrix::beProductTOf(const DoubleMatrix& a, const DoubleMatrix& b, double scale)
{
	DoubleMatrix tmp(b.Tanspose());
	*this = a * tmp;
	this->Scale(scale);
}


void DoubleMatrix::plusProductSymm(const DoubleMatrix& a, const DoubleMatrix& b, double scale)
{
	if (this->isEmpty())
	{
		this->nRows = a.nColumns;
		this->nColumns = b.nColumns;
		this->values.assign(a.nColumns * b.nColumns, 0.);
	}
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = i; j <= nColumns; j++)
		{
			double sum = 0;
			for (int k = 1; k <= a.nRows; k++)
			{
				sum += a(k, i) * b(k, j);
			}
			(*this)(i, j) += sum * scale;
		}
	}

}



void DoubleMatrix::plusProductUnSymm(const DoubleMatrix& a, const DoubleMatrix& b, double scale)
{
	if (this->isEmpty())
	{
		this->nRows = a.nColumns;
		this->nColumns = b.nColumns;
		this->values.assign(a.nColumns * b.nColumns, 0.);
	}
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			double sum = 0;
			for (int k = 1; k <= a.nRows; k++)
			{
				sum += a(k, i) * b(k, j);
			}
			(*this)(i, j) += sum * scale;
		}
	}
}

void DoubleMatrix::symmetrized()
{
	for (int i = 2; i <= nRows; i++)
	{
		for (int j = 1; j < i; j++)
		{
			(*this)(i, j) = (*this)(j, i);
		}
	}
}

void DoubleMatrix::beUnitMatrix()
{
	if (this->nColumns != this->nRows)
	{
		ERROR("can not be unit matrix!");
	}
	this->zero();
	for (int i = 1; i <= this->nRows; i++)
	{
		(*this)(i, i) = 1.;
	}
}

void DoubleMatrix::beMatrixForm(const DoubleArray& src)
{
	
	if (src.GetSize() == 9)
	{
		this->resize(3, 3);
		(*this)(1, 1) = src(1);
		(*this)(2, 2) = src(2);
		(*this)(3, 3) = src(3);
		(*this)(2, 3) = src(4);
		(*this)(1, 3) = src(5);
		(*this)(1, 2) = src(6);
		(*this)(3, 2) = src(7);
		(*this)(3, 1) = src(8);
		(*this)(2, 1) = src(9);
	}
	else if (src.GetSize() == 6)
	{
		this->resize(3, 3);
		(*this)(1, 1) = src(1);
		(*this)(2, 2) = src(2);
		(*this)(3, 3) = src(3);
		(*this)(2, 3) = src(4);
		(*this)(1, 3) = src(5);
		(*this)(1, 2) = src(6);
		(*this)(3, 2) = src(4);
		(*this)(3, 1) = src(5);
		(*this)(2, 1) = src(6);
	}
	else if (src.GetSize() == 4)
	{
		this->resize(2, 2);
		(*this)(1, 1) = src(1);
		(*this)(2, 2) = src(2);
		(*this)(1, 2) = src(3);
		(*this)(2, 1) = src(4);
	}
	else if (src.GetSize() == 3)
	{
		this->resize(2, 2);
		(*this)(1, 1) = src(1);
		(*this)(2, 2) = src(2);
		(*this)(1, 2) = src(3);
		(*this)(2, 1) = src(3);
	}
	else
	{
		ERROR("unknown vector size!");
	}
}

void DoubleMatrix::beSubMatrix(const DoubleMatrix& src, const IntArray& a, const IntArray& b)
{
	this->resize(a.getSize(), b.getSize());
	for (int i = 1; i <= a.getSize(); i++)
	{
		for (int j = 1; j <= b.getSize(); j++)
		{
			(*this)(i, j) = src(a(i), b(j));
		}
	}
}

DoubleMatrix::DoubleMatrix(std::initializer_list<std::initializer_list<double>>mat)
{
	this->resize(mat.size(), mat.begin()->size());
	auto p = this->values.begin();
	for (auto col : mat)
	{
		for (auto x : col)
		{
			*p = x;
			p++;
		}
	}
}
DoubleMatrix& DoubleMatrix::operator=(std::initializer_list<std::initializer_list<double>> mat)
{
	this->resize(mat.size(), mat.begin()->size());
	auto p = this->values.begin();
	for (auto col : mat)
	{
		for (auto x : col)
		{
			*p = x;
			p++;
		}
	}
	return *this;
}

void DoubleMatrix::printfYourself() const
{
	printf("(%d x %d): \n", nRows, nColumns);
	if (nRows <= 250 && nColumns <= 250) {
		for (int i = 1; i <= nRows; ++i) {
			for (int j = 1; j <= nColumns && j <= 100; ++j) {
				printf("%10.3e  ", (*this)(i, j));
			}

			printf("\n");
		}
	}
	else {
		for (int i = 1; i <= nRows && i <= 20; ++i) {
			for (int j = 1; j <= nColumns && j <= 10; ++j) {
				printf("%10.3e  ", (*this)(i, j));
			}
			if (nColumns > 10) printf(" ...");
			printf("\n");
		}
		if (nRows > 20)  printf(" ...\n");
	}
}

void DoubleMatrix::copySubVectoRow(const DoubleArray& src, int r, int c)
{
	c--;
	int cols = src.GetSize();
	int nr = r;
	int nc = c + cols;
	if (this->GetRow() < nr || this->GetColumn() < nc)
	{
		this->resizewithDatas(max_fem(this->GetRow(), nr), max_fem(this->GetColumn(), nc));
	}

	for (int i = 1; i <= cols; i++)
	{
		(*this)(nr, c + i) = src(i);
	}
}

void DoubleMatrix::resizewithDatas(int r, int c)
{// 旧数据会被保存的
	// 如果是原矩阵大小就不能加入
	if (r == this->nRows&&c == this->nColumns)
	{
		return;
	}

	// move 这个命令比较好用
	DoubleMatrix old(*this);
	this->nRows = r;
	this->nColumns = c;
	this->values.resize(r * c);

	int ii = min_fem(r, old.GetRow());
	int jj = min_fem(c, old.GetColumn());
	for (int i = 1; i <= ii; i++)
	{
		for (int j = 1; j <= jj; j++)
		{
			(*this)(i, j) = old(i, j);
		}
	}

}

void DoubleMatrix::rotatedwith(const DoubleMatrix& r, char mode)
{
	DoubleMatrix tmp;
	if (mode == 'n')
	{
		tmp.TProduct(r, *this);
		this->beProductOf(tmp, r);
	}
	else if (mode == 't')
	{
		tmp.beProductOf(r, *this);
		this->beProductTOf(tmp, r);
	}
	else
	{
		ERROR("unsupported mode!");
	}
}

void DoubleMatrix::setSubMatrix(const DoubleMatrix& src, int r, int c)
{
	r--;
	c--;
	int srcrow = src.GetRow();
	int srccol = src.GetColumn();
	for (int i = 1; i <= srcrow; i++)
	{
		for (int j = 1; j <= srccol; j++)
		{
			(*this)(r + i, c + j) = src(i, j);
		}
	}
}

void DoubleMatrix::assemble(const DoubleMatrix& src, const IntArray& loc)
{
	int r, c, size, ii, jj;
	r = src.GetRow();
	c = src.GetColumn();
	size = loc.getSize();
	if (r != c || r != size || c != size)
	{
		ERROR("cannot assemble! size not match");
	}
	for (int i = 1; i <= r; i++)
	{
		if ((ii = loc(i)))
		{
			for (int j = 1; j <= c; j++)
			{
				if ((jj = loc(j)))
				{
					(*this)(ii, jj) += src(i, j);
				}
			}
		}
		
	}
}

bool DoubleMatrix::solveForRhs(const DoubleArray& b, DoubleArray& answer, bool transpose)
{//@to debug
	if (!this->isSquare())
	{
		ERROR("cannot solve");
	}
	if (nRows != b.GetSize())
	{
		ERROR("dimension mismatch!");
	}
	int pivRow;
	double piv, linkomb, help;
	DoubleMatrix trans;
	DoubleMatrix *K;
	if (transpose)
	{
		trans = this->Tanspose();
		K = &trans;
	}
	else
	{
		K = this;
	}
	answer = b;
	for (int i = 1; i < nRows ; ++i)
	{
		piv = fabs((*K)(i, i));
		pivRow = i;
		for (int j = i + 1; j <= nRows; ++j)
		{
			if (fabs((*K)(j, i)) > piv)
			{
				pivRow = j;
				piv = fabs((*K)(j, i));
			}
		}
		if (piv < 1.0e-20)
		{
			return false;
		}
		if (pivRow != i)
		{
			for (int j = i; j <= nRows; j++)
			{
				help = (*K)(i, j);
				(*K)(i, j) = (*K)(pivRow, j);
				(*K)(pivRow, j) = help;
			}
			help = answer(i);
			answer(i) = answer(pivRow);
			answer(pivRow) = help;
		}
		for (int j = i + 1; j <= nRows; j++)
		{
			linkomb = (*K)(j, i) / (*K)(i, i);
			for (int k = i; k <= nRows; k++)
			{
				(*K)(j, k) -= (*K)(i, k) * linkomb;
			}
			answer(j) -= answer(i) * linkomb;
		}
	}
	for (int i = nRows; i >= 1; i--)
	{
		help = 0.;
		for (int j = i + 1; j <= nRows; j++)
		{
			help += (*K)(i, j) * answer(j);
		}
		answer(i) = (answer(i) - help) / (*K)(i, i);
	}
	return true;
}

void DoubleMatrix::beDiagonal(const DoubleArray& src)
{
	int n = src.GetSize();
	this->resize(n, n);
	this->zero();
	for (int i = 1; i <= n; i++)
	{
		(*this)(i, i) = src(i);
	}
}

void DoubleMatrix::be_I_X_I_Matrix()
{
	this->resize(6, 6);
	this->zero();
	values[0] = values[1] = values[2] = 1.;
	values[6] = values[7] = values[8] = 1.;
	values[12] = values[13] = values[14] = 1.;
}

void DoubleMatrix::beIs()
{
	this->resize(6, 6);
	this->zero();
	values[0] = values[7] = values[14] = 1.;
	values[21] = values[28] = values[35] = 0.5;
}

void DoubleMatrix::beId()
{
	this->resize(6, 6);
	this->zero();
	values[0] = values[7] = values[14] = 2. / 3.;
	values[1] = values[2] = values[6] = values[8] = values[12] = values[13] = -1. / 3.;
	values[21] = values[28] = values[35] = 0.5;
}

void DoubleMatrix::beDyadicProductOf(const DoubleArray& a, const DoubleArray& b)
{
	this->resize(a.GetSize(), b.GetSize());
	for (int i = 1; i <= a.GetSize(); i++)
	{
		for (int j = 1; j <= b.GetSize(); j++)
		{
			(*this)(i, j) = a(i) * b(j);
		}
	}

}

void DoubleMatrix::plusDyadicSymm(const DoubleArray& a, double scale)
{
	if (this->isEmpty())
	{
		this->nRows = a.GetSize();
		this->nColumns = a.GetSize();
		this->values.assign(this->nRows * this->nColumns, 0.);
	}
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = i; j <= nColumns; j++)
		{
			(*this)(i, j) += a(i) * a(j) * scale;
		}
	}
}

void DoubleMatrix::plusDyadicUnsymm(const DoubleArray& a, double scale)
{
	if (this->isEmpty())
	{
		this->nRows = a.GetSize();
		this->nColumns = a.GetSize();
		this->values.assign(this->nRows * this->nColumns, 0.);
	}
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			(*this)(i, j) += a(i) * a(j) * scale;
		}
	}
}

void DoubleMatrix::beSymPartOf(const DoubleMatrix& src)
{
	if (src.nRows != src.nColumns)
	{
		ERROR("src is not square!");
	}
	this->resize(src.nRows, src.nColumns);
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			(*this)(i, j) = 0.5 * (src(i, j) + src(j, i));
		}
	}
}

void DoubleMatrix::beAntiSymPartOf(const DoubleMatrix& src)
{
	if (src.nRows != src.nColumns)
	{
		ERROR("src is not square!");
	}
	this->resize(src.nRows, src.nColumns);
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			(*this)(i, j) = 0.5 * (src(i, j) - src(j, i));
		}
	}
}

void DoubleMatrix::eigProjectionOf(const DoubleMatrix& vec)
{
	if (vec.GetRow() == 2)
	{// 11 22 12
		this->resize(3, 2);
		(*this)(1, 1) = vec(1, 1) * vec(1, 1);
		(*this)(1, 2) = vec(1, 2) * vec(1, 2);

		(*this)(2, 1) = vec(2, 1) * vec(2, 1);
		(*this)(2, 2) = vec(2, 2) * vec(2, 2);

		(*this)(3, 1) = vec(1, 1) * vec(2, 1);
		(*this)(3, 2) = vec(1, 2) * vec(2, 2);
	}
	else if (vec.GetRow() == 3)
	{ // 11 22 33 23 13 12
		this->resize(6, 3);
		(*this)(1, 1) = vec(1, 1) * vec(1, 1);
		(*this)(1, 2) = vec(1, 2) * vec(1, 2);
		(*this)(1, 3) = vec(1, 3) * vec(1, 3);

		(*this)(2, 1) = vec(2, 1) * vec(2, 1);
		(*this)(2, 2) = vec(2, 2) * vec(2, 2);
		(*this)(2, 3) = vec(2, 3) * vec(2, 3);

		(*this)(3, 1) = vec(3, 1) * vec(3, 1);
		(*this)(3, 2) = vec(3, 2) * vec(3, 2);
		(*this)(3, 3) = vec(3, 3) * vec(3, 3);

		(*this)(3, 1) = vec(3, 1) * vec(3, 1);
		(*this)(3, 2) = vec(3, 2) * vec(3, 2);
		(*this)(3, 3) = vec(3, 3) * vec(3, 3);

		(*this)(4, 1) = vec(2, 1) * vec(3, 1);
		(*this)(4, 2) = vec(2, 2) * vec(3, 2);
		(*this)(4, 3) = vec(2, 3) * vec(3, 3);

		(*this)(5, 1) = vec(1, 1) * vec(3, 1);
		(*this)(5, 2) = vec(1, 2) * vec(3, 2);
		(*this)(5, 3) = vec(1, 3) * vec(3, 3);

		(*this)(6, 1) = vec(1, 1) * vec(2, 1);
		(*this)(6, 2) = vec(1, 2) * vec(2, 2);
		(*this)(6, 3) = vec(1, 3) * vec(2, 3);
	}
}

DoubleMatrix DoubleMatrix::cutoff(DoubleMatrix A, int i, int j)
{
	DoubleMatrix B(A.GetRow() - 1, A.GetColumn() - 1);
	for (int c = 1; c <= B.GetRow(); c++)
	{
		for (int r = 1; r <= B.GetColumn(); r++)
		{
			B(c, r) = A(c + (c >= i), r + (r >= j));
		}
	}
	return B;
}

double DoubleMatrix::det(DoubleMatrix A)
{
	if (A.GetColumn() == 1)
	{
		return A(1, 1);
	}
	
	else
	{
		double ans = 0;
		for (int j = 1; j <= A.GetColumn(); j++)
		{
			ans += A(1, j) *det(cutoff(A, 1, j)) * (j % 2 ? -1 : 1);
		}
		return ans;
	}
}

DoubleMatrix DoubleMatrix::company(DoubleMatrix A)
{
	DoubleMatrix B(A.GetRow(), A.GetColumn());
	for (int i = 1; i <= B.GetRow(); i++)
	{
		for (int j = 1; j <= B.GetColumn(); j++)
		{
			B(j, i) = det(cutoff(A, i, j)) * ((i + j) % 2 ? -1 : 1);
		}
	}
	return B;
}

DoubleMatrix DoubleMatrix::num_mul(DoubleMatrix A, double num)
{
	DoubleMatrix B(A.GetRow(), A.GetColumn());
	for (int i = 1; i <= B.GetRow(); i++)
	{
		for (int j = 1; j <= B.GetColumn(); j++)
		{
			B(i, j) = A(i, j) * num;
		}
	}
	return B;
}

void DoubleMatrix::matrix_inverse(DoubleMatrix& a, DoubleMatrix& b)
{
	int i, j, k;
	int n = a.GetRow();
	double max, temp;
	DoubleMatrix t(n, n);//定义一个临时矩阵
	//将a矩阵临时存放在矩阵t中
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
		{
			t(i, j) = a(i, j);

		}
	}
	//初始化B矩阵为单位矩阵
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
		{
			b(i, j) = (i == j) ? (double)1 : 0;

		}
	}
	//进行列主消元，找到每一列的主元；
	for (i = 1; i <= n; i++)
	{
		max = t(i, i);
		k = i;
		for (j = i + 1; j <= n; j++)
		{
			if (fabs(t(j, i)) > fabs(max))
			{
				max = t(j, i);
				k = j;
			}

		}
		//cout<<"the max number is" <<max<<endl;
		//如果主元所在行不是第i行，则进行 行交换
		if (k != i)
		{
			for (j = 1; j <= n; j++)
			{
				temp = t(i, j);
				t(i, j) = t(k, j);
				t(k, j) = temp;
			}

		}
		if (t(i, i) == 0)
		{
			cout << "\nthe matrix does not exist inverse matrix\n";
			break;
		}
		//获取列主元素
		temp = t(i, i);

		//将主元所在的行进行单元化处理
		//cout << "\nthe temp is "<<temp<<endl;
		for (j = 1; j <= n; j++)
		{
			t(i, j) = t(i, j) / temp;;
			b(i, j) = b(i, j) / temp;
		}
		for (j = 1; j <= n; j++)
		{
			if (j != i)
			{
				temp = t(j, i);
				for (k = 1; k <= n; k++)
				{
					t(j, k) = t(j, k) - temp * t(i, k);
					b(j, k) = b(j, k) - temp * b(i, k);
				}
			}
		}


	}

}


void DoubleMatrix::PlusProduct(const DoubleMatrix& a, const DoubleMatrix& b, double scale)
{
	if (this->isEmpty())
	{
		this->nRows = a.GetRow();
		this->nColumns = b.GetColumn();
		this->values.assign(a.GetRow() * b.GetColumn(), 0);
	}
	for (int i = 1; i <= nRows; i++)
	{
		for (int j = 1; j <= nColumns; j++)
		{
			double sum = 0.;
			for (int k = 1; k <= a.GetColumn(); k++)
			{
				sum += a(i, k) * b(k, j);
			}
			(*this)(i, j) += sum * scale;
		}
	}
}

void DoubleMatrix::copyMtrxOf(const DoubleMatrix& src, int rstart, int cstart, int rend, int cend)
{
	int nr = rend - rstart + 1;
	int nc = cend - cstart + 1;
	this->resize(nr, nc);
	this->zero();
	for (int i = 1; i <= nr; i++)
	{
		for (int j = 1; j <= nc; j++)
		{
			(*this)(i, j) = src(rstart + i - 1, cstart + j - 1);
		}
	}
}

/* void DoubleMatrix::beTanspose()
{

	for (int i = 1; i <= this->nRows; i++)
	{
		for (int j = i + 1; j <= this->nColumns; j++)
		{
			double tmp;
			tmp = (*this)(i, j);
			(*this)(i, j) = (*this)(j, i);
			(*this)(j, i) = tmp;
		}
	}
}*/

//DoubleMatrix cutoff(DoubleMatrix A, int i, int j)
//{
//	DoubleMatrix B(A.GetRow() - 1, A.GetColumn() - 1);
//	for (int c = 1; c <= B.GetRow(); c++)
//	{
//		for (int r = 1; r <= B.GetColumn(); r++)
//		{
//			B(c, r) = A(c + (c >= i), r + (r >= j));
//		}
//	}
//	return B;
//}
//
//double det(DoubleMatrix A)
//{
//	if (A.GetColumn() == 1)
//	{
//		return A(1, 1);
//	}
//	double ans = 0;
//	for (int j = 1; j <= A.GetColumn(); j++)
//	{
//		ans += A(1, j) * det(cutoff(A, 1, j)) * (j % 2 ? -1 : 1);
//	}
//	return ans;
//}
//
//DoubleMatrix company(DoubleMatrix A)
//{
//	DoubleMatrix B(A.GetRow(), A.GetColumn());
//	for (int i = 1; i <= B.GetRow(); i++)
//	{
//		for (int j = 1; j <= B.GetColumn(); j++)
//		{
//			B(j, i) = det(cutoff(A, i, j)) * ((i + j) % 2 ? -1 : 1);
//		}
//	}
//	return B;
//}
//
//DoubleMatrix num_mul(DoubleMatrix A, double num)
//{
//	DoubleMatrix B(A.GetRow(), A.GetColumn());
//	for (int i = 1; i <= B.GetRow(); i++)
//	{
//		for (int j = 1; j <= B.GetColumn(); j++)
//		{
//			B(i, j) = A(i, j) * num;
//		}
//	}
//	return B;
//}

DoubleMatrix DoubleMatrix::inv44(DoubleMatrix t)
{
	int i, j, k;
	DoubleMatrix s(4, 4);
	for (i = 1; i <= 3; i++)
	{
		int pivot = i;
		double pivotsize = t(i, i);
		if (pivotsize < 0)
		{
			pivotsize = -pivotsize;
		}
		for (j = i + 1; j <= 4; j++)
		{
			double temp = t(j, i);
			if (temp < 0) { temp = -temp; }
			if (temp > pivotsize) { pivot = j; pivotsize = temp; }
		}
		if (pivotsize == 0) { return DoubleMatrix(); }
		if (pivot != i)
		{
			for (j = 1; j <= 4; j++)
			{
				double tmp;
				tmp = t(i, j);
				t(i, j) = t(pivot , j);
				t.SetValues(pivot , j, tmp);
				tmp = s(i, j);
				s(i, j) = s(pivot , j);
				s(pivot, j) = tmp;
			}
		}
		for (j = i + 1; j <= 4; j++)
		{
			double f = t(j, i) / t(i, i);
			for (k = 1; k <= 4; k++)
			{
				t(j, k) -= t(i, k) * f;
				s(j, k) -= s(i, k) * f;
			}
		}
	}
	for (i = 4; i >= 1; i--)
	{
		double f;
		f = t(i, i);
		if (f == 0) {
			return DoubleMatrix();
		}
		for (j = 1; j <= 4; j++)
		{
			t(i, j) /= f;
			s(i, j) /= f;
		}
		for (j = 1; j < i; j++)
		{
			f = t(j, i);
			for (k = 1; k <= 4; k++)
			{
				t(j, k) -= f * t(i, k);
				s(j, k) -= f * s(i, k);
			}
		}
	}
	return s;
}
DoubleMatrix DoubleMatrix::inv(DoubleMatrix& a)
{
	if (a.GetRow() != a.GetColumn())
	{
		ERROR("求逆矩阵不为方阵");
		return NULL;
	}
	int* is, * js, i, j, k;
	int n = a.GetRow();
	double temp, fmax;
	DoubleMatrix p(a.GetRow(), a.GetColumn());
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n; j++)
		{
			p(i, j) = a(i, j);
		}
	}
	is = new int[n];
	js = new int[n];
	for (k = 0; k < n; k++)
	{
		fmax = 0.0;
		for (i = k; i < n; i++)
			for (j = k; j < n; j++)
			{
				temp = std::fabs(p(i + 1, j + 1));
				if (temp > fmax)
				{
					fmax = temp;
					is[k] = i; js[k] = j;
				}
			}
		if ((fmax + 1.0) == 1.0)
		{
			delete[]is;
			delete[]js;
			ERROR("不存在逆矩阵");
			return NULL;
		}
		if ((i = is[k]) != k)
		{
			for (j = 0; j < n; j++)
			{
				swap2(p(k + 1, j + 1), p(i + 1, j + 1));
			}
		}
		if ((j = js[k]) != k)
		{
			for (i = 0; i < n; i++)
			{
				swap2(p(i + 1, k + 1), p(i + 1, j + 1));
			}
		}
		p(k + 1, k + 1) = 1.0 / p(k + 1, k + 1);
		for (j = 0; j < n; j++)
		{
			if (j != k)
			{
				p(k + 1, j + 1) *= p(k + 1, k + 1);
			}
		}
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				for (j = 0; j < n; j++)
				{
					if (j != k)
					{
						p(i + 1, j + 1) = p(i + 1, j + 1) - p(i + 1, k + 1) * p(k + 1, j + 1);
					}
				}
			}
		}
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				p(i + 1, k + 1) *= -p(k + 1, k + 1);
			}
		}
	}
	
	for (k = n - 1; k >= 0; k--)
	{
		if ((j = js[k]) != k)
		{
			for (i = 0; i < n; i++)
			{
				swap2(p(j + 1, i + 1), p(k + 1, i + 1));
			}
		}
		if ((i = is[k]) != k)
		{
			for (j = 0; j < n; j++)
			{
				swap2(p(j + 1, i + 1), p(j + 1, k + 1));
			}
		}
	}
	delete[]is;
	delete[]js;
	*this= p;
	return *this;
}
void DoubleMatrix::swap2(double& a, double& b)
{
	double c;
	c = a;
	a = b;
	b = c;
}
