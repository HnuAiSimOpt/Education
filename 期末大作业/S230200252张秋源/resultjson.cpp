#include "tetra4-struct.hpp"
#include<nlohmann/json.hpp>
#include<vector>
#include<array>
using namespace std;
using json = nlohmann::json;

/*
Put calculated data to json file for further analysis.
*/
json outputTetra4json(vector<Point3DResult>& p3dRstVec, vector<Tetra4Result>& t4RstVec){
    json nodeContent;
    int nodecount = 0;
    for(vector<Point3DResult>::iterator it = p3dRstVec.begin(); it != p3dRstVec.end(); it++){
        json nd;
        nd["nodenum"] = nodecount;
        nd["position"] = {it->x, it->y, it->z};
        nd["displace"] = {it->u, it->v, it->w};
        nodeContent.push_back(nd);
        nodecount++;
    }
    json elementContent;
    for(vector<Tetra4Result>::iterator it = t4RstVec.begin(); it != t4RstVec.end(); it++){
        json el;
        el["node"] = it->node;
        el["number"] = it->elementNumber;
        el["strain"] = it->strainVec;
        el["stress"] = it->stressVec;
        el["vonmises"] = it->vonMises;
        elementContent.push_back(el);
    }
    json out;
    out["nodes"] = nodeContent;
    out["elements"] = elementContent;
    return out;
}
/*
To generate FEM analysis result as json output like:
start node number must be 0
{
    "nodes":[{
        "position":[x0, y0, z0],
        "displace":[u0, v0, w0]
    },
    {
        "position":[x1, y1, z1],
        "displace":[u1, v1, w1]
    }],
    "elements":[{
        "node":[n1, n2, n3, n4],
        "number":x,
        "strain":[s1, s2, ..., s6],
        "stress":[s1, s2, ..., s6]
    },
    {
        ...
    }
    ]
}
*/
