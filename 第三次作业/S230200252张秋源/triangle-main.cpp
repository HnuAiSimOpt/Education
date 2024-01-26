#include "triangle.cpp"

int main(int argc, char** argv){
    /*
    Generate a 12 * 5 rectangle, divide to 231 points, make it 21 * 11,
    each element size is 0.6 * 0.5 triangle
    */
   int pointNum = 231;
   vector<Point2D> elementPoint;
   elementPoint.reserve(pointNum);
   for(int i = 0; i < 21; i++){
    for(int j = 0; j < 11; j++){
        elementPoint.push_back(Point2D(0.6 * i, 0.5 * j));
    }
   }
   /*
   Define triangle elements
   */
  vector<TrianglePoint> triLabel;
  for(int i = 0; i < 20; i++){
    for(int j = 0; j < 10; j++){
        triLabel.push_back(TrianglePoint(i * 11 + j, (i + 1) * 11 + j, (i + 1) * 11 + j + 1));
        triLabel.push_back(TrianglePoint(i * 11 + j, (i + 1) * 11 + j + 1, i * 11 + j + 1));
    }
  }
    /*
    Set young's modulus to 1e2
    Set poisson ratio to 0.2
    */
    Eigen::Matrix3d dmat = getDMat(1e2, 0.2);

    /*
    Set the displacement of the first 11 points on the most left to zero
    */
    vector<DisplaceBound> dispBound;
    for(int i = 0; i < 11; i++){
        dispBound.push_back(DisplaceBound(i, 0, 0));
    }
    /*
    Apply force to node 220 to 230
    */
    vector<double> force(2 * pointNum, 0);
    for(int i = 0; i < 11; i++){
        force[(220 + i) * 2] =  0.1 * i; //x direction point force
        force[(220 + i) * 2 + 1] = 0.1 * i; //y direction point force
    }
    /*
    Get the attribute for each element
    */
    vector<TriangleElement> triEle = getTriangleElement(elementPoint, triLabel, dmat);
    /*
    Assemble K matrix
    */
    Eigen::SparseMatrix<double, Eigen::RowMajor> kmat = assembleKMat(triEle, pointNum);
    /*
    Apply displacement constrain
    */
    applyDisplacementConstrain(kmat, force, dispBound);
    vector<double> displacements = solveDisplacement(kmat, force);
    vector<PointResult> pointRst = genPointResult(elementPoint, displacements);
    cout << "Final solution:\n";
    printPointResult(cout, pointRst);
    cout.flush();
    return 0;

}

