#include<iostream>
#include<vector>
#include<array>
#include<exception>
#include<stdexcept>
#include<cmath>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/SuperLUSupport>
#include"tetra4-struct.hpp" //Where structures are defined
using namespace std;
/*
Eigen library vector to std vector
*/
vector<double> eigenToStdVec(Eigen::VectorXd& evec){ //Eigen vector to std vector
    vector<double> out(evec.data(), evec.data() + evec.size());
    return out;
}
/*
Get D Matrix for a specified material
*/
Eigen::Matrix<double, 6, 6> getDMat(Material& m){ // Get D matrix
    Eigen::Matrix<double, 6, 6> outmat = Eigen::Matrix<double, 6, 6>::Zero();
    double a = m.modulusE / (1.0 + m.poissonRatio) / (1.0 - 2 * m.poissonRatio);
    outmat(0, 0) = 1.0 - m.poissonRatio;
    outmat(1, 1) = 1.0 - m.poissonRatio; 
    outmat(2, 2) = 1.0 - m.poissonRatio; 
    outmat(3, 3) = (1.0 - 2 * m.poissonRatio) / 2;
    outmat(4, 4) = (1.0 - 2 * m.poissonRatio) / 2;
    outmat(5, 5) = (1.0 - 2 * m.poissonRatio) / 2;
    outmat(0, 1) = m.poissonRatio;
    outmat(0, 2) = m.poissonRatio;
    outmat(1, 0) = m.poissonRatio;
    outmat(1, 2) = m.poissonRatio;
    outmat(2, 0) = m.poissonRatio;
    outmat(2, 1) = m.poissonRatio;
    outmat *= a;
    return outmat;
}
/*
double TetraVolume(array<double, 4>& x, array<double, 4>& y, array<double, 4>& z){
    Eigen::Matrix<double, 4, 4> m = Eigen::Matrix<double, 4, 4>::Ones();
    for(int i = 0; i < 4; i++){
        m(i, 1) = x[i];
        m(i, 2) = y[i];
        m(i, 3) = z[i];
    }
    return m.determinant() / 6;
}
*/


/*
Get B matrix and K matrix for each element
*/
vector<Tetra4Element> getBKmats(vector<Point3D>& pt, vector<Tetra4NoBK>& ele, Material& materi){
    Eigen::Matrix<double, 6, 6> getDMat(Material&);

    Eigen::Matrix<double, 6, 6> dmat = getDMat(materi); //Get D matrix from material
    array<double, 4> x; //x position of 4 points
    array<double, 4> y; //y position of 4 points
    array<double, 4> z; //z position of 4 points
    Eigen::Matrix<double, 4, 4> m = Eigen::Matrix<double, 4, 4>::Ones(); //Matrix for calculating the volume
    Eigen::Matrix<double, 4, 4> invm; //Inverse of the matrix above
    array<double, 3> bcd; //Storing the value for generating B matrix
    Eigen::Matrix<double, 6, 12> bmat = Eigen::Matrix<double, 6, 12>::Zero(); //B Matrix
    Eigen::Matrix<double, 12, 12> kmat;
    vector<Tetra4Element> out;
    double vol; //Volume of the tetrahydron
    for(vector<Tetra4NoBK>::iterator it = ele.begin(); it != ele.end(); it++){ 
        for(int i = 0; i < 4; i++){ //Read position of 4 points
            x[i] = pt[it->node[i]].x;
            y[i] = pt[it->node[i]].y;
            z[i] = pt[it->node[i]].z;
            m(i, 1) = x[i];
            m(i, 2) = y[i];
            m(i, 3) = z[i];
        }
        vol = m.determinant() / 6;
        invm = m.inverse();
        for(int i = 0; i < 4; i++){ //fill B matrix
            bcd[0] = invm(1, i);
            bcd[1] = invm(2, i);
            bcd[2] = invm(3, i);
            bmat(0, 3 * i) = bcd[0];
            bmat(1, 3 * i + 1) = bcd[1];
            bmat(2, 3 * i + 2) = bcd[2];
            bmat(3, 3 * i) = bcd[1];
            bmat(3, 3 * i + 1) = bcd[0];
            bmat(4, 3 * i + 1) = bcd[2];
            bmat(4, 3 * i + 2) = bcd[1];
            bmat(5, 3 * i) = bcd[2];
            bmat(5, 3 * i + 2) = bcd[0];
        }
        kmat = bmat.transpose() * dmat * bmat * vol;
        Tetra4Element outele(*it, bmat, kmat);
        out.push_back(outele);
        //it->bmat = bmat;
        //it->kmat = bmat.transpose() * dmat * bmat * vol; //Calculate K matrix
    }
    return out;
}

/*
Assemble these K matrix
*/
Eigen::SparseMatrix<double, Eigen::RowMajor> assembleKmat(vector<Tetra4Element>& ele, int pointNum){
    Eigen::SparseMatrix<double, Eigen::RowMajor> mat(pointNum * 3, pointNum * 3);
    vector<Eigen::Triplet<double>> tplist;
    for(vector<Tetra4Element>::iterator it = ele.begin(); it != ele.end(); it++){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                tplist.push_back(Eigen::Triplet<double>(it->node[i] * 3, it->node[j] * 3));
                tplist.push_back(Eigen::Triplet<double>(it->node[i] * 3 + 1, it->node[j] * 3));
                tplist.push_back(Eigen::Triplet<double>(it->node[i] * 3 + 2, it->node[j] * 3));
                tplist.push_back(Eigen::Triplet<double>(it->node[i] * 3, it->node[j] * 3 + 1));
                tplist.push_back(Eigen::Triplet<double>(it->node[i] * 3 + 1, it->node[j] * 3 + 1));
                tplist.push_back(Eigen::Triplet<double>(it->node[i] * 3 + 2, it->node[j] * 3 + 1));
                tplist.push_back(Eigen::Triplet<double>(it->node[i] * 3, it->node[j] * 3 + 2));
                tplist.push_back(Eigen::Triplet<double>(it->node[i] * 3 + 1, it->node[j] * 3 + 2));
                tplist.push_back(Eigen::Triplet<double>(it->node[i] * 3 + 2, it->node[j] * 3 + 2));
            }
        }
    }
    mat.setFromTriplets(tplist.begin(), tplist.end());
    for(vector<Tetra4Element>::iterator it = ele.begin(); it != ele.end(); it++){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                mat.coeffRef(it->node[i] * 3, it->node[j] * 3) += it->kmat(3 * i, 3 * j);
                mat.coeffRef(it->node[i] * 3 + 1, it->node[j] * 3) += it->kmat(3 * i + 1, 3 * j);
                mat.coeffRef(it->node[i] * 3 + 2, it->node[j] * 3) += it->kmat(3 * i + 2, 3 * j);
                mat.coeffRef(it->node[i] * 3, it->node[j] * 3 + 1) += it->kmat(3 * i, 3 * j + 1);
                mat.coeffRef(it->node[i] * 3 + 1, it->node[j] * 3 + 1) += it->kmat(3 * i + 1, 3 * j + 1);
                mat.coeffRef(it->node[i] * 3 + 2, it->node[j] * 3 + 1) += it->kmat(3 * i + 2, 3 * j + 1);
                mat.coeffRef(it->node[i] * 3, it->node[j] * 3 + 2) += it->kmat(3 * i, 3 * j + 2);
                mat.coeffRef(it->node[i] * 3 + 1, it->node[j] * 3 + 2) += it->kmat(3 * i + 1, 3 * j + 2);
                mat.coeffRef(it->node[i] * 3 + 2, it->node[j] * 3 + 2) += it->kmat(3 * i + 2, 3 * j + 2);
            }
        }
    }
    return mat;
}
/*
Use set-1 method to apply displacement constrain
*/
void applyDisplacementConstrain(Eigen::SparseMatrix<double, Eigen::RowMajor>& akmat, vector<double>& forceVec, vector<DisplaceBound>& dispBound){
    for(vector<DisplaceBound>::iterator it = dispBound.begin(); it != dispBound.end(); it++){
        /*
        Use "Set 1 method"
        */
        akmat.row(it->p * 3) *= 0;
        akmat.row(it->p * 3 + 1) *= 0;
        akmat.row(it->p * 3 + 2) *= 0;
        akmat.coeffRef(it->p * 3, it->p * 3) = 1;
        akmat.coeffRef(it->p * 3 + 1, it->p * 3 + 1) = 1;
        akmat.coeffRef(it->p * 3 + 2, it->p * 3 + 2) = 1;
        forceVec[it->p * 3] = it->u;
        forceVec[it->p * 3 + 1] = it->v;
        forceVec[it->p * 3 + 2] = it->w;
    }
}

/*
Solve displacement for each point, Use another library superLU
*/
vector<double> solveDisplacement(Eigen::SparseMatrix<double, Eigen::RowMajor>& akmat, vector<double>& forceVec){ 
    vector<double> eigenToStdVec(Eigen::VectorXd&);

    Eigen::VectorXd out;
    Eigen::SuperLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
    Eigen::VectorXd forceMat = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(forceVec.data(), forceVec.size());
    solver.compute(akmat);
    if(solver.info() != Eigen::Success){
        cerr << "Solving failed" << endl;
        out = Eigen::VectorXd::Zero(1); //failure only return 1 value
    }
    else{
        out = solver.solve(forceMat);
    }
    vector<double> out1 = eigenToStdVec(out);
    return out1;
}
/*
Put the result of displacement vector to the struct of each point
*/
vector<Point3DResult> getPointResult(vector<Point3D>& p3dv, vector<double>& displaceVec){
    if(p3dv.size() * 3 != displaceVec.size()){
        throw runtime_error("displaceVec wrong size");
    }
    vector<Point3DResult> out;
    out.reserve(p3dv.size());
    size_t ptSize = p3dv.size();
    for(size_t i = 0; i < ptSize; i++){
        out.push_back(Point3DResult(p3dv[i], displaceVec[i * 3], displaceVec[i * 3 + 1], displaceVec[i * 3 + 2]));
    }
    return out;
}
/*
Get the result of tetra4 elements
*/
vector<Tetra4Result> getTetra4Result(vector<Point3DResult>& pointResult, vector<Tetra4Element>& ele, Material& materi){
    Eigen::Matrix<double, 12, 1> uMat;
    Eigen::Matrix<double, 6, 1> epsilonMat;
    Eigen::Matrix<double, 6, 1> sigmaMat;
    Eigen::Matrix<double, 6, 6> dmat = getDMat(materi); //Get D matrix from material
    vector<Tetra4Result> out;
    Tetra4Result t4rst;
    for(vector<Tetra4Element>::iterator it = ele.begin(); it != ele.end(); it++){
        for(int i = 0; i < 4; i++){ //Generate displacement matrix U
            uMat[3 * i] = pointResult[it->node[i]].u;
            uMat[3 * i + 1] = pointResult[it->node[i]].v;
            uMat[3 * i + 2] = pointResult[it->node[i]].w;
        }
        epsilonMat = it->bmat * uMat; //Calculate strain
        sigmaMat = dmat * epsilonMat; //Calculate stress
        t4rst.node = it->node;
        t4rst.elementNumber = it->elementNumber;
        for(int i = 0; i < 6; i++){
            t4rst.strainVec[i] = epsilonMat[i];
            t4rst.stressVec[i] = sigmaMat[i];
        }
        Eigen::Matrix<double, 3, 3> stressTensor;
        stressTensor(0, 0) = t4rst.stressVec[0];
        stressTensor(1, 1) = t4rst.stressVec[1];
        stressTensor(2, 2) = t4rst.stressVec[2];
        stressTensor(0, 1) = t4rst.stressVec[3];
        stressTensor(1, 0) = t4rst.stressVec[3];
        stressTensor(0, 2) = t4rst.stressVec[4];
        stressTensor(2, 0) = t4rst.stressVec[4];
        stressTensor(1, 2) = t4rst.stressVec[5];
        stressTensor(2, 1) = t4rst.stressVec[5];
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> solver(stressTensor);
        if(solver.info() != Eigen::Success) throw runtime_error("Stress tensor matrix solving error");
        Eigen::Vector3d pStress = solver.eigenvalues(); //principle stress
        t4rst.vonMises = sqrt(0.5 * ((pStress[0] - pStress[1]) * (pStress[0] - pStress[1]) + 
        (pStress[1] - pStress[2]) * (pStress[1] - pStress[2]) + (pStress[0] - pStress[2])
        * (pStress[0] - pStress[2])));
        out.push_back(t4rst);
    }
    return out;
}