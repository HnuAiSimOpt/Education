#include "LinearStatic.h"
#include <iostream>
#include "DoubleMatrix.h"
using namespace std;
#include "DoubleArray.h"
#include<string>
#include<fstream>
#include<stdio.h>


/**
 * \brief 学号：B2302S0105	姓名：刘文广
 * \param 线性三角形
 */

void  fematiso(DoubleMatrix& p)
{
	double elastic = 1.0;
	double poisson = 0.0;
	DoubleMatrix u(3, 3); u.zero();
	u.SetValues(1, 1, 1);		u.SetValues(1, 2, poisson);
	u.SetValues(2, 1, poisson);	u.SetValues(2, 2, 1);
	u.SetValues(3, 3, (1 - poisson) / 2);

	for (int i = 1; i <= 3; i++)
	{
		for (int j = 1; j <= 3; j++)
		{
			p(i, j) = (elastic / (1 - poisson * poisson)) * u(i, j);

		}
	}

}

IntArray feeldof(IntArray& p, int a, int b)
{
	int edof = a * b;
	int k = 0;
	IntArray index(6);
	for (int i = 1; i <= a; ++i)
	{
		int start = (p.at(i - 1) - 1) * b;
		for (int j = 1; j <= b; ++j)
		{
			k = k + 1;
			index.SetValue(k, (start + j));

		}
	}
	return index;
}

DoubleMatrix fekine2d(int nnel, DoubleArray& p1, DoubleArray& p2)
{
	DoubleMatrix km2(3, 6);
	for (int i = 1; i <= nnel; i++)
	{
		int i1 = (i - 1) * 2 + 1;
		int i2 = i1 + 1;

		km2.SetValues(1, i1, p1.at(i - 1));
		km2.SetValues(2, i2, p2.at(i - 1));
		km2.SetValues(3, i1, p2.at(i - 1));
		km2.SetValues(3, i2, p1.at(i - 1));

	}
	return km2;
}



void feasmb_2(DoubleMatrix& p1, DoubleMatrix& p2, IntArray& p3)
{
	int edof = p3.getSize();
	int ii = 0, jj = 0;
	for (int i = 1; i <=edof; i++)
	{
		ii = p3.at(i - 1);
		for (int j = 1; j <= edof; ++j)
		{
			jj = p3.at(j - 1);
			p1.SetValues(ii, jj, p1(ii, jj) + p2(i, j));
		}

	}
}

void feaplyc2(DoubleMatrix& p1, DoubleArray& p2, IntArray& p3, IntArray& p4)
{
	int n = p3.getSize();
	int sdof = p1.GetColumn();
	for (int i = 1; i <= n; i++)
	{
		int c = p3(i);
		for (int j = 1; j <= sdof; j++)
		{
			p1.SetValues(c, j, 0);
		}
		p1.SetValues(c, c, 1.0);
		double d = p4(i);
		p2.SetValue(c, d);

	}


}
int main(int argc, char* argv[])
{
	Domain* d = new Domain;
	LinearStatic ls(d);
	ls.assemble();

	//input data for control parameters

	int lengthx = 4;              //length of x - axis side of problem
	int lengthy = 2;             // length of y - axis side of problem

	int lx = 16;                   // number of element in x - axis
	int ly = 8;                    // number of element in y - axis

	int nel = 2 * lx * ly;          // number of element
	int nnel = 3;                   // number of nodes per element

	int ndof = 2;                 // number of dofs per node
	int nnode = (lx + 1) * (ly + 1);//  total number of nodes in system
	int sdof = nnode * ndof;         // total system dofs
	int edof = nnel * ndof;          // degrees of freedom per element

	int emodule = 1;					// elastic modulus
	int poisson = 0;					// Poisson's ratio

	int fload = -1;						// the total load


	//input data for nodal coordinate values
	//gcoord(i, j) where i->node no.and j->x or y

	DoubleMatrix x0(153, 2); x0.zero();
	int num = 0;
	for (int i = 1; i <= lx + 1; i++)
	{
		for (int j = 1; j <= ly + 1; j++)
		{
			++num;
			x0.SetValues(num, 1, (double)(i - 1) * (double)lengthx / (double)lx);
			x0.SetValues(num, 2, (double)((0 - 0.5) * (double)lengthy * (double)(1 + (double)(lx + 1 - i) / (double)lx) * (double)((1 - (double)(j - 1) / (double)ly))));

		}
	}
	/*input data for nodal connectivity for each element
	 nodes(i, j) where i->element no.and j->connected nodes*/
	num = 0;
	vector<vector<int>>nodes(256, vector<int>(0));
	for (int i = 1; i <= lx; ++i)
	{
		for (int j = 1; j <= ly; ++j)
		{
			nodes[num].push_back((ly + 1) * (i - 1) + j);
			nodes[num].push_back((ly + 1) * i + j);
			nodes[num].push_back((ly + 1) * (i - 1) + j + 1);
			++num;
			nodes[num].push_back((ly + 1) * i + j);
			nodes[num].push_back((ly + 1) * i + j + 1);
			nodes[num].push_back((ly + 1) * (i - 1) + j + 1);
			++num;

		}
	}
	/*DoubleMatrix nodes(256, 3); nodes.zero();
	for (int i = 1; i <= lx; ++i)
	{
		for (int j = 1; j <= ly; ++j)
		{
			++num;
			nodes.SetValues(num, 1, (ly + 1) * (i - 1) + j);
			nodes.SetValues(num, 2, (ly + 1) * i + j);
			nodes.SetValues(num, 3, (ly + 1) * (i - 1) + j + 1);

			++num;
			nodes.SetValues(num, 1, (ly + 1) * i + j);
			nodes.SetValues(num, 2, (ly + 1) * i + j + 1);
			nodes.SetValues(num, 3, (ly + 1) * (i - 1) + j + 1);
		}
	}*/



	//input data for boundary conditions

	IntArray bcdof(0);
	IntArray bcval(0);
	for (int i = 1; i <= ly + 1; i++)
	{
		bcdof.push_back(1 + 2 * (i - 1));
		bcdof.push_back(2 + 2 * (i - 1));
		bcval.push_back(0);
		bcval.push_back(0);

	}

	//initialization of matricesand vectors  创建一系列矩阵

	DoubleArray ff(sdof);	ff.zero();			//system force vector


	DoubleMatrix k(edof, edof);				//initialization of element matrix
	k.zero();


	DoubleMatrix kk(sdof, sdof);				//system matrix
	kk.zero();


	DoubleArray disp(sdof); disp.zero();					//system displacement vector


	DoubleArray eldisp(edof);				//element displacement vector
	eldisp.zero();


	DoubleMatrix stress(nel, 3);				//matrix containing stress components
	stress.zero();

	DoubleMatrix strain(nel, 3);				//matrix containing stress components
	strain.zero();

	IntArray index(edof);				//index vector
	index.zero();

	vector<vector<int>> kinmtx(3, vector<int>(edof, 0));				//kinematic matrix

	DoubleMatrix matmtx(3, 3);				//constitutive matrix
	matmtx.zero();

	/*------------------------------------------------ -
		 compute element matrices and vectors and assemble
	 ------------------------------------------------ -*/

	fematiso(matmtx);

	for (int iel = 1; iel <= nel; iel++)
	{
		IntArray nd(3); nd.zero();

		int  nd1 = nodes[iel - 1][0];
		int  nd2 = nodes[iel - 1][1];
		int  nd3 = nodes[iel - 1][2];

		nd.SetValue(1, nd1);
		nd.SetValue(2, nd2);
		nd.SetValue(3, nd3);

		double x1, y1, x2, y2, x3, y3;
		x1 = x0.GetValues(nd1, 1);      y1 = x0.GetValues(nd1, 2);
		x2 = x0.GetValues(nd2, 1);      y2 = x0.GetValues(nd2, 2);
		x3 = x0.GetValues(nd3, 1);      y3 = x0.GetValues(nd3, 2);

		index = feeldof(nd, nnel, ndof);

		//y1 = fabs(y1); y2 = fabs(y2); y3 = fabs(y3);

		double area = 0.5 * (x1 * y2 + x2 * y3 + x3 * y1 - x1 * y3 - x2 * y1 - x3 * y2);
		double area2 = area * 2;

		DoubleArray dhdx(3); dhdx.zero();
		DoubleArray dhdy(3); dhdy.zero();

		dhdx.SetValue(1, (1 / area2) * (y2 - y3));
		dhdx.SetValue(2, (1 / area2) * (y3 - y1));
		dhdx.SetValue(3, (1 / area2) * (y1 - y2));

		dhdy.SetValue(1, (1 / area2) * (x3 - x2));
		dhdy.SetValue(2, (1 / area2) * (x1 - x3));
		dhdy.SetValue(3, (1 / area2) * (x2 - x1));



		//vector<vector<int>> kinmtx2(3, vector<int>(6, 0));
		//kinmtx2 = fekine2d(nnel, dhdx, dhdy);
		//vector<vector<int>> kinmtx2T(6, vector<int>(3, 0));
		//for (int i = 1; i <= 3; i++)
		//{
		//	for (int j = 1; j <= 6; j++)
		//	{
		//		kinmtx2T[j][i] = kinmtx2[i][j];
		//	}
		//}
		//k = kinmtx2T * matmtx * kinmtx2 * area;

		DoubleMatrix kinmtx2(3, 6); kinmtx2.zero();
		kinmtx2 = fekine2d(nnel, dhdx, dhdy);
		DoubleMatrix kinmtx2T(6, 3); kinmtx2T.zero();
		for (int i = 1; i <= 3; i++)
		{
			for (int j = 1; j <= 6; j++)
			{
				kinmtx2T.SetValues(j, i, kinmtx2.GetValues(i, j));
			}
		}
		k = kinmtx2T * matmtx * kinmtx2 * area;      //element stiffness matrice

		feasmb_2(kk, k, index);

	}

	DoubleMatrix kk1 = kk;
	//创建负载矩阵，每个节点的负载为-1
	 //ff.assign(-1.0);
	ff.SetValue(306, -1.0);
	feaplyc2(kk, ff, bcdof, bcval);	// 对总刚度矩阵和负载矩阵进行约束化 即左边第一列单元位移为0

							//kk.solveForRhs(ff,disp);    // disp = kk\ff;	求出各个节点的位移
	DoubleMatrix b(306, 306); b.zero();
	kk.matrix_inverse(kk, b);
	/*int n = b.GetColumn();
	*/
	double temp = 0;
	for (int i = 1; i <= b.GetRow(); i++)
	{
		temp = 0;
		for (int j = 1; j <= b.GetColumn(); j++)
		{
			temp += b(i, j) * ff(j);
		}
		disp.SetValue(i, temp);
	}
	
	double EU = 0.5 * disp.calenergy(kk1);
	
	DoubleArray U_fem = disp;

	double energy = 0;
	double energy0 = 0;
	double denergy = 0;

	for (int ielp = 1; ielp <= nel; ielp++)
	{
		IntArray nd(3); nd.zero();

		int  nd1 = nodes[ielp - 1][0];
		int  nd2 = nodes[ielp - 1][1];
		int  nd3 = nodes[ielp - 1][2];

		nd.SetValue(1, nd1);
		nd.SetValue(2, nd2);
		nd.SetValue(3, nd3);

		double x1, y1, x2, y2, x3, y3;
		x1 = x0.GetValues(nd1, 1);      y1 = x0.GetValues(nd1, 2);
		x2 = x0.GetValues(nd2, 1);      y2 = x0.GetValues(nd2, 2);
		x3 = x0.GetValues(nd3, 1);      y3 = x0.GetValues(nd3, 2);

		// 单元中心(重心）
		double xcentre = (x1 + x2 + x3) / 3;
		double ycentre = (y1 + y2 + y3) / 3;

		index = feeldof(nd, nnel, ndof);

		for (int i = 1; i <= edof; i++)
		{

			double temp = disp.at(index.at(i - 1) - 1);
			eldisp.SetValue(i, temp);

		}
		double area = 0.5 * (x1 * y2 + x2 * y3 + x3 * y1 - x1 * y3 - x2 * y1 - x3 * y2);
		double area2 = area * 2;

		DoubleArray dhdx(3); dhdx.zero();
		DoubleArray dhdy(3); dhdy.zero();

		dhdx.SetValue(1, (1 / area2) * (y2 - y3));
		dhdx.SetValue(2, (1 / area2) * (y3 - y1));
		dhdx.SetValue(3, (1 / area2) * (y1 - y2));

		dhdy.SetValue(1, (1 / area2) * (x3 - x2));
		dhdy.SetValue(2, (1 / area2) * (x1 - x3));
		dhdy.SetValue(3, (1 / area2) * (x2 - x1));

		DoubleMatrix kinmtx2(3, 6); kinmtx2.zero();
		kinmtx2 = fekine2d(nnel, dhdx, dhdy);

		DoubleArray estrain(3); estrain.zero();
		DoubleArray estress(3); estrain.zero();

		estrain.MtrxDotVect(kinmtx2, eldisp);
		estress.MtrxDotVect(matmtx, estrain);

		for (int i = 1; i <= 3; i++)
		{
			strain(ielp, i) = estrain(i);
			stress(ielp, i) = estress(i);
		}

		energy += 0.5 * estrain.calenergy(matmtx) * area;

	}

	vector<vector<int>>neigh_node; neigh_node.resize(153);
	IntArray indneigh(153); indneigh.zero();
	for (int i = 1; i <= nel; i++)
	{
		for (int j = 1; j <= 3; j++)
		{
			++indneigh(nodes[i - 1][j - 1]);
			neigh_node[nodes[i - 1][j - 1] - 1].push_back(i);
		}
	}


	DoubleMatrix stress_node(3, nnode); stress_node.zero();

	for (int inode = 1; inode <= nnode; inode++)
	{
		int numel = indneigh(inode);
		for (int i = 1; i <= numel; i++)
		{
			int ind_nel = neigh_node[inode - 1].at(i - 1);
			for (int j = 1; j <= 3; j++)
			{
				stress_node(j, inode) = stress_node(j, inode) + stress(ind_nel, j);
			}

		}
		stress_node(1, inode) = stress_node(1, inode) / numel;
		stress_node(2, inode) = stress_node(2, inode) / numel;
		stress_node(3, inode) = stress_node(3, inode) / numel;
	}

		FILE* fid_out2;
		fid_out2 = fopen("result_beam01.plt", "w");
		fprintf(fid_out2, "TITLE=\"test case governed by poisson equation\"\n");
		fprintf(fid_out2, "VARIABLES=\"x\" \"y\" \"u\" \"v\" \"sigax\"  \"sigmay\" \"sigmaxy\"\n");
		fprintf(fid_out2, "ZONE T=\"flow-field\", N= %8d,E=%8d,ET=TRIANGLE, F=FEPOINT\n", nnode, nel);
		for (int i = 1; i <= nnode; i++)
		{
			fprintf(fid_out2, " %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e\n", x0(i, 1), x0(i, 2), disp(2 * i - 1), disp(2 * i), stress_node(1, i), stress_node(2, i), stress_node(3, i));
		}

		for (int i = 1; i <= nel; i++)
		{
			fprintf(fid_out2, "%8d %8d %8d\n", nodes[i - 1][0], nodes[i - 1][1], nodes[i - 1][2]);
		}
		fclose(fid_out2);


}
