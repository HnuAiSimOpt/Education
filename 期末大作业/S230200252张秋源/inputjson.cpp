#include "tetra4-struct.hpp"
#include<nlohmann/json.hpp>
#include<vector>
#include<array>
#include<exception>
#include<stdexcept>
using namespace std;
using json = nlohmann::json;

vector<Point3D> getPoints(json fj){
    vector<Point3D> out;
    if(!fj.contains("nodes")) throw runtime_error("Json file not containing node");
    json nodesContent = fj["nodes"];
    if(!nodesContent.is_array()) throw runtime_error("Nodes type error");
    for(json::iterator it = nodesContent.begin(); it != nodesContent.end(); it++){
        Point3D p3d;
        json k = *it;
        if(!(k.is_array()) || !(k.size() == 3)) throw runtime_error("Not array in nodes");
        p3d.x = k[0].get<double>();
        p3d.y = k[1].get<double>();
        p3d.z = k[2].get<double>();
        out.push_back(p3d);
    }
    return out;
}

vector<Tetra4NoBK> getElements(json fj){
    vector<Tetra4NoBK> out;
    if(!fj.contains("elements")) throw runtime_error("Json file not containing elements");
    json elementContent = fj["elements"];
    if(!elementContent.is_array()) throw runtime_error("elements attribute is not array");
    for(json::iterator it = elementContent.begin(); it != elementContent.end(); it++){
        Tetra4NoBK nbk;
        json k = *it;
        if(!(k.is_object()) || !(k.contains("number") || !(k.contains("node")))) throw runtime_error("element attribute not valid");
        nbk.elementNumber = k["number"].get<int>();
        json nd = k["node"];
        if(!(nd.is_array())) throw runtime_error("element attribute not valid");
        for(int i = 0; i < 4; i++){
            nbk.node[i] = nd[i].get<int>();
        }
        out.push_back(nbk);
    }
    return out;
}

vector<DisplaceBound> getConstrains(json fj){
    vector<DisplaceBound> out;
    if(!fj.contains("constrains")) throw runtime_error("Json file not containing constrains");
    json dbd = fj["constrains"];
    if(!dbd.is_array()) throw runtime_error("Constrain type error");
    for(json::iterator it = dbd.begin(); it != dbd.end(); it++){
        DisplaceBound d;
        json k = *it;
        if(!(k.contains("node")) || !(k.contains("displace"))) throw runtime_error("Json file constrain type error");
        d.p = k["node"].get<int>();
        json dsp = k["displace"];
        if(!(dsp.is_array())) throw runtime_error("Json file constrain type error");
        d.u = dsp[0].get<double>();
        d.v = dsp[1].get<double>();
        d.w = dsp[2].get<double>();
        out.push_back(d);
    }
    return out;
}

vector<double> getForceVec(json fj, int pointNum){
    vector<double> forceVec(pointNum * 3, 0.0);
    if(!fj.contains("forces")) throw runtime_error("Json file not containing forces");
    json forceContent = fj["forces"];
    if(!forceContent.is_array()) throw runtime_error("Json file force type error");
    for(json::iterator it = forceContent.begin(); it != forceContent.end(); it++){
        json k = *it;
        if(!(k.contains("node")) || !(k.contains("forceval"))) throw runtime_error("Json file force type error");
        int n = k["node"].get<int>();
        json fcval = k["forceval"];
        if(!(fcval.is_array())) throw runtime_error("Json file force type error");
        forceVec[3 * n] += fcval[0].get<double>();
        forceVec[3 * n + 1] += fcval[1].get<double>();
        forceVec[3 * n + 2] += fcval[2].get<double>();
    }
    return forceVec;
}

Material getMaterial(json fj){
    if(!fj.contains("material")) throw runtime_error("Json file not containing material");
    json materialContent = fj["material"];
    Material m;
    if(!(materialContent.contains("modulus-of-elasticity")) || !(materialContent.contains("poisson-ratio"))) throw runtime_error("Material format not correct");
    m.modulusE = materialContent["modulus-of-elasticity"].get<double>();
    m.poissonRatio = materialContent["poisson-ratio"].get<double>();
    return m;
}
/*
Read from a file
That file must contains position of the nodes, elements, displacement constrains, forces.
{
    "nodes":[
        [x0, y0, z0],
        [x1, y1, z1],
        ...
    ],
    "elements":[
        {"number":5000, "node":[n1, n2, n3, n4]},
        {"number":5001, "node":[n1, n2, n3, n4]},
        ...
    ],
   "constrains":[
        {"node":10, "displace":[x1, y1, z1]},
        {"node":11, "displace":[x2, y2, z2]},
        ...
   ],
   "forces":[
        {"node":100, "forceval":[1.23, 4.56, 7.89]},
        {"node":101, "forceval":[1.23, 4.56, 7.89]},
        ...
   ],
   "material":{
    "modulus-of-elasticity":1e3,
    "poisson-ratio":0.3
   }
}
*/