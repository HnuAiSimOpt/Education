#ifndef INPUTJSON_HPP
#define INPUTJSON_HPP

#include "tetra4-struct.hpp"
#include<nlohmann/json.hpp>
#include<vector>
#include<array>
#include<exception>
#include<stdexcept>
using namespace std;
using json = nlohmann::json;

vector<Point3D> getPoints(json fj);

vector<Tetra4NoBK> getElements(json fj);

vector<DisplaceBound> getConstrains(json fj);

vector<double> getForceVec(json fj, int pointNum);

Material getMaterial(json fj);
#endif