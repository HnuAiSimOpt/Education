#include<iostream>
#include<nlohmann/json.hpp>
#include<vector>
#include<map>
#include<string>
#include<fstream>
#include<chrono>
#include"inputjson.hpp"
#include"resultjson.hpp"
#include"tetra4.hpp"
#include"argparse.hpp"

using namespace std;
using json = nlohmann::json;

int main(int argc, char** argv){
    map<string, vector<string>> command = argparse(argc, argv);
    if(command.find("-i") == command.end()) throw runtime_error("Not setting input tag -i");
    if(command["-i"].size() == 0) throw runtime_error("Input file name not specified");
    string inputfilename = command["-i"][0];
    if(command.find("-o") == command.end()) throw runtime_error("Not setting input tag -o");
    if(command["-o"].size() == 0) throw runtime_error("Output file name not specified");
    string outputfilename = command["-o"][0];
    ifstream inputfile(inputfilename.data());
    if(!inputfile.is_open()) throw runtime_error("Input file not able to open");
    json fj = json::parse(inputfile);

    vector<Point3D> pt = getPoints(fj);
    int pointNum = pt.size();
    vector<Tetra4NoBK> t4nbk = getElements(fj);
    vector<DisplaceBound> dispBound = getConstrains(fj);
    vector<double> forceVec = getForceVec(fj, pointNum);
    Material mt = getMaterial(fj);

    auto startime = chrono::high_resolution_clock::now();

    vector<Tetra4Element> t4Ele = getBKmats(pt, t4nbk, mt);
    Eigen::SparseMatrix<double, Eigen::RowMajor> akmat = assembleKmat(t4Ele, pointNum);
    applyDisplacementConstrain(akmat, forceVec, dispBound);
    vector<double> displaceVec = solveDisplacement(akmat, forceVec);
    vector<Point3DResult> p3drst = getPointResult(pt, displaceVec);
    vector<Tetra4Result> t4rst = getTetra4Result(p3drst, t4Ele, mt);

    auto stoptime = chrono::high_resolution_clock::now();
    cout << "Calculation finished" << endl;
    cout << "Nodes: " << p3drst.size() << endl;
    cout << "Elements: " << t4rst.size() << endl;
    cout << "Time used: " << chrono::duration_cast<chrono::milliseconds>(stoptime - startime) << endl;

    json outjson = outputTetra4json(p3drst, t4rst);
    ofstream outputfile(outputfilename.data());
    if(!outputfile.is_open()) throw runtime_error("Output file can't be opened");
    outputfile << outjson << endl;
    inputfile.close();
    outputfile.close();
    return 0;
}