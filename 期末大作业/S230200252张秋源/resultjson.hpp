#ifndef RESULTJSON_HPP
#define RESULTJSON_HPP

#include "tetra4-struct.hpp"
#include<nlohmann/json.hpp>
using namespace std;
using json = nlohmann::json;

/*
Put calculated data to json file for further analysis.
*/
json outputTetra4json(vector<Point3DResult>& p3dRstVec, vector<Tetra4Result>& t4RstVec);
#endif