#ifndef TETRA4_HPP
#define TETRA4_HPP
#include<iostream>
#include<vector>
#include<array>
#include<exception>
#include<stdexcept>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include"tetra4-struct.hpp"
using namespace std;

/*
Eigen library vector to std vector
*/
vector<double> eigenToStdVec(Eigen::VectorXd& evec);

/*
Get D Matrix for a specified material
*/
Eigen::Matrix<double, 6, 6> getDMat(Material& m);

/*
Get B matrix and K matrix for each element
*/
vector<Tetra4Element> getBKmats(vector<Point3D>& pt, vector<Tetra4NoBK>& ele, Material& materi);

/*
Assemble these K matrix
*/
Eigen::SparseMatrix<double, Eigen::RowMajor> assembleKmat(vector<Tetra4Element>& ele, int pointNum);

/*
Use set-1 method to apply displacement constrain
*/
void applyDisplacementConstrain(Eigen::SparseMatrix<double, Eigen::RowMajor>& akmat, vector<double>& forceVec, vector<DisplaceBound>& dispBound);

/*
Solve displacement for each point, Use another library superLU
*/
vector<double> solveDisplacement(Eigen::SparseMatrix<double, Eigen::RowMajor>& akmat, vector<double>& forceVec);

/*
Put the result of displacement vector to the struct of each point
*/
vector<Point3DResult> getPointResult(vector<Point3D>& p3dv, vector<double>& displaceVec);

/*
Get the result of tetra4 elements
*/
vector<Tetra4Result> getTetra4Result(vector<Point3DResult>& pointResult, vector<Tetra4Element>& ele, Material& materi);

#endif