/*
Program made by Qiuyuan Zhang, from Hunan University
S230200252
*/

#include<iostream>
#include<vector>
#include<string>
#include<Eigen/Dense>
#include<Eigen/Sparse>

using namespace std;
/*
A point
*/
class Point2D{
    public:
        double x;
        double y;
        Point2D(double, double);
        Point2D();
};
Point2D::Point2D(double xval, double yval){
    x = xval;
    y = yval;
}
Point2D::Point2D(){
    x = 0;
    y = 0;
}

/*
Point, combined with its displacement result
*/
class PointResult : public Point2D{
    public:
        double u;
        double v;
};
/*
The labels of triangle points
*/
class TrianglePoint{
    public:
        int p[3];
        TrianglePoint(int, int, int);
        TrianglePoint();
};
TrianglePoint::TrianglePoint(int p0, int p1, int p2){
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
}
TrianglePoint::TrianglePoint(){
    p[0] = 0;
    p[1] = 0;
    p[2] = 0;
}
/*
Inherit from trianglePoint, contains area,
B matrix and K matrix
*/
class TriangleElement : public TrianglePoint{
    public:
        double area;
        Eigen::Matrix<double, 3, 6> bmat;
        Eigen::Matrix<double, 6, 6> kmat;
};
/*
Displacement constrain 
*/
class DisplaceBound{
    public:
        int p;
        double u;
        double v;

        DisplaceBound(int, double, double);
};

DisplaceBound::DisplaceBound(int point, double xdisplace, double ydisplace) : p(point), u(xdisplace), v(ydisplace){

}
/*
Get 3x3 Eigen Matrix D by young's modulus E and
poisson ratio mu
*/
vector<double> eigenToStdVec(Eigen::VectorXd& evec){
    vector<double> out(evec.data(), evec.data() + evec.size());
    return out;
}
Eigen::Matrix3d getDMat(double e, double mu){
    double a =  e / (1.0 + mu) / (1.0 - 2 * mu);
    Eigen::Matrix3d outmat = Eigen::Matrix3d::Zero();
    outmat(0, 0) = 1.0 - mu;
    outmat(1, 1) = 1.0 - mu;
    outmat(0, 1) = mu;
    outmat(1, 0) = mu;
    outmat(2, 2) = (1.0 - 2 * mu) / 2;
    outmat *= a;
    return outmat;
}

/*
Return a vector of area of those units
*/
vector<double> getArea(vector<Point2D>& pxy, vector<TrianglePoint>& ps){
    vector<double> out;
    out.reserve(ps.size());
    Eigen::Matrix3d a = Eigen::Matrix3d::Zero();
    double b;
    for(vector<TrianglePoint>::iterator i = ps.begin(); i != ps.end(); i++){
        a(0, 0) = 1;
        a(1, 0) = 1;
        a(2, 0) = 1;
        a(0, 1) = pxy[i->p[0]].x;
        a(1, 1) = pxy[i->p[1]].x;
        a(2, 1) = pxy[i->p[2]].x;
        a(2, 0) = pxy[i->p[0]].y;
        a(2, 1) = pxy[i->p[1]].y;
        a(2, 2) = pxy[i->p[2]].y;
        b = a.determinant() / 2;
        out.push_back(b);
    }
    return out;
}

/*
Return a vector of 3 by 6 matrix of B matrix
*/
vector<Eigen::Matrix<double, 3, 6>> getBMat(vector<Point2D>& pxy, vector<TrianglePoint>& ps, vector<double>& area){
    vector<Eigen::Matrix<double, 3, 6>> out;
    out.reserve(ps.size());
    Eigen::Matrix<double, 3, 6> m;
    double b1, c1, b2, c2, b3, c3;
    int counter1 = 0;
    for(vector<TrianglePoint>::iterator i = ps.begin(); i != ps.end(); i++){
        m = Eigen::Matrix<double, 3, 6>::Zero();
        b1 = pxy[i->p[1]].y - pxy[i->p[2]].y;
        c1 = pxy[i->p[2]].x - pxy[i->p[1]].x;
        b2 = pxy[i->p[2]].y - pxy[i->p[0]].y;
        c2 = pxy[i->p[0]].x - pxy[i->p[2]].x;
        b3 = pxy[i->p[0]].y - pxy[i->p[1]].y;
        c3 = pxy[i->p[1]].x - pxy[i->p[0]].x;
        m(0, 0) = b1;
        m(1, 1) = c1;
        m(2, 0) = c1;
        m(2, 1) = b1;
        m(0, 2) = b2;
        m(1, 3) = c2;
        m(2, 2) = c2;
        m(2, 3) = b2;
        m(0, 4) = b3;
        m(1, 5) = c3;
        m(2, 4) = c3;
        m(2, 5) = b3;
        m *= 1.0 / area[counter1] / 2;
        counter1++;
        out.push_back(m);
    }
    return out;
}

/*
Get the status of each element,
Return a vector of TriangleElement class
*/
vector<TriangleElement> getTriangleElement(vector<Point2D>& pxy, vector<TrianglePoint>& ps, Eigen::Matrix3d& dmat){
    vector<double> getArea(vector<Point2D>&, vector<TrianglePoint>&);
    vector<Eigen::Matrix<double, 3, 6>> getBMat(vector<Point2D>&, vector<TrianglePoint>&, vector<double>&);
    
    vector<double> areas = getArea(pxy, ps);
    vector<Eigen::Matrix<double, 3, 6>> bmats = getBMat(pxy, ps, areas);
    vector<TriangleElement> out;
    out.reserve(ps.size());
    TriangleElement ele;
    for(int i = 0; i < ps.size(); i++){
        ele.p[0] = ps[i].p[0];
        ele.p[1] = ps[i].p[1];
        ele.p[2] = ps[i].p[2];
        ele.area = areas[i];
        ele.bmat = bmats[i];
        ele.kmat = areas[i] * bmats[i].transpose() * dmat * bmats[i];
        out.push_back(ele);
    }
    return out;

}

/*
Assemble SparseMatrix
*/
Eigen::SparseMatrix<double, Eigen::RowMajor> assembleKMat(vector<TriangleElement>& triEle, int pointNum){
    Eigen::SparseMatrix<double, Eigen::RowMajor> mat(pointNum * 2, pointNum * 2);
    vector<Eigen::Triplet<double>> tplist;
    for(vector<TriangleElement>::iterator i = triEle.begin(); i != triEle.end(); i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                tplist.push_back(Eigen::Triplet<double>(i->p[j] * 2, i->p[k] * 2, i->kmat(j * 2, k * 2)));
                tplist.push_back(Eigen::Triplet<double>(i->p[j] * 2 + 1, i->p[k] * 2, i->kmat(j * 2 + 1, k * 2)));
                tplist.push_back(Eigen::Triplet<double>(i->p[j] * 2, i->p[k] * 2 + 1, i->kmat(j * 2, k * 2 + 1)));
                tplist.push_back(Eigen::Triplet<double>(i->p[j] * 2 + 1, i->p[k] * 2 + 1, i->kmat(j * 2 + 1, k * 2 + 1)));
            }
        }
    }
    mat.setFromTriplets(tplist.begin(), tplist.end());
    return mat;
}

/*
Apply displacement constrain by a class called DisplaceBound
Use "set 1 method"
*/
void applyDisplacementConstrain(Eigen::SparseMatrix<double, Eigen::RowMajor>& akmat, vector<double>& forcemat, vector<DisplaceBound>& disBound){
    int row = akmat.rows();
    int col = akmat.cols();
    int targetpt;
    for(vector<DisplaceBound>::iterator i = disBound.begin(); i != disBound.end(); i++){
        targetpt = i->p;
        akmat.row(2 * targetpt) *= 0;
        akmat.coeffRef(2 * targetpt, 2 * targetpt) = 1;
        forcemat[2 * targetpt] = i->u;
        akmat.row(2 * targetpt + 1) *= 0;
        akmat.coeffRef(2 * targetpt + 1, 2 * targetpt + 1) = 1;
        forcemat[2 * targetpt + 1] = i->v;
    }
}

/*
Solve the displacements of all points as vector
*/
vector<double> solveDisplacement(Eigen::SparseMatrix<double, Eigen::RowMajor>& akmat, vector<double>& force){
    vector<double> eigenToStdVec(Eigen::VectorXd&);
    
    Eigen::VectorXd out;
    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>> solver;
    Eigen::VectorXd forcemat = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(force.data(), force.size());
    solver.compute(akmat);
    if(solver.info() != Eigen::Success){
        cerr << "Solving failed" << endl;
        out = Eigen::VectorXd::Zero(1); //failure only return 1 value
    }
    else{
        out = solver.solve(forcemat);
    }
    vector<double> out1 = eigenToStdVec(out);
    return out1;
}

/*
collect the calculating result of each point
*/
vector<PointResult> genPointResult(vector<Point2D>& p2d, vector<double>& displace){
    vector<PointResult> out;
    PointResult r1;
    out.reserve(p2d.size());
    for(int i = 0; i < p2d.size(); i++){
        r1.x = p2d[i].x;
        r1.y = p2d[i].y;
        r1.u = displace[2 * i];
        r1.v = displace[2 * i + 1];
        out.push_back(r1);
    }
    return out;
}

/*
put the result to ostream buffer
*/
void printPointResult(ostream& ost, vector<PointResult>& pointRst){
    int counter1 = 0;
    for(vector<PointResult>::iterator i = pointRst.begin(); i != pointRst.end(); i++){
        ost << "Point label:" << counter1 << '\n';
        ost << "x=" << i->x << ' ' << "y=" << i->y << ' ';
        ost << "u=" << i->u << ' ' << "v=" << i->v << '\n';
        counter1++;
    }
}