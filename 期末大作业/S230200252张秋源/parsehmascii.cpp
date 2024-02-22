#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<nlohmann/json.hpp>
#include"argparse.hpp"
using namespace std;
using json = nlohmann::json;

//Not finished yet

vector<string> fileToStrings(const string& filename){
    vector<string> outstrs;
    ifstream f1;
    string tmpstring;
    f1.open(filename, ios::in);
    if(f1.is_open()){
        f1.seekg(0, f1.end);
        size_t fsize = f1.tellg();
        f1.seekg(0, f1.beg);
        while (f1.tellg() < fsize)
        {
            getline(f1, tmpstring);
            outstrs.push_back(tmpstring);
        }
        f1.close();
    }
    else{
        cerr << "The file " << filename << "can't be opened" << endl;
    }
    return outstrs;
}
/*
Separate a string by comma, within bracket
*/
vector<string> sepString(const string& ist, char lbrack, char rbrack, char sep){
    vector<string> out;
    size_t p1 = ist.find(lbrack);
    size_t p2 = ist.find(rbrack);
    if(p1 != string::npos && p2 != string::npos && p2 > p1){
        string s2 = ist.substr(p1 + 1, p2 - p1 - 1);
        size_t currPos = 0;
        size_t nextPos = s2.find(sep);
        while(true){
            if(nextPos == string::npos){
                out.push_back(s2.substr(currPos, s2.size() - currPos));
                break;
            }
            out.push_back(s2.substr(currPos, nextPos - currPos));
            currPos = nextPos + 1;
            nextPos = s2.find(sep, currPos);
        }
    }
    return out;
}

json parseHMAscii(string& filename){
    vector<string> vst = fileToStrings(filename);
    json allContent;
    json nodeContent;
    json tetra4Content;
    json constrainContent;
    json forceContent;
    for(vector<string>::iterator it = vst.begin(); it != vst.end(); it++){
        if(it->starts_with("*node")){
            json a3;
            vector<string> sp = sepString(*it, '(', ')', ',');
            a3 = {stod(sp[1]), stod(sp[2]), stod(sp[3])};
            nodeContent.push_back(a3);
        }
        if(it->starts_with("*constraint")){
            json a4;
            vector<string> sp = sepString(*it, '(', ')', ',');
            a4["node"] = (stoi(sp[2]) - 1);
            a4["displace"] = {stod(sp[4]), stod(sp[5]), stod(sp[6])};
            constrainContent.push_back(a4);
        }
        if(it->starts_with("*tetra4")){
            json a4;
            vector<string> sp = sepString(*it, '(', ')', ',');
            a4["number"] = (stoi(sp[0]));
            a4["node"] = {stoi(sp[2]) - 1, stoi(sp[3]) - 1, stoi(sp[4]) - 1, stoi(sp[5]) - 1};
            tetra4Content.push_back(a4);
        }
        if(it->starts_with("*force")){
            json a4;
            vector<string> sp = sepString(*it, '(', ')', ',');
            a4["node"] = (stoi(sp[2]) - 1);
            a4["forceval"] = {stod(sp[4]), stod(sp[5]), stod(sp[6])};
            forceContent.push_back(a4);
        }
    }
    allContent["nodes"] = nodeContent;
    allContent["elements"] = tetra4Content;
    allContent["constrains"] = constrainContent;
    allContent["forces"] = forceContent;
    return allContent;
}

int main(int argc, char** argv){
    map<string, vector<string>> command = argparse(argc, argv);
    if(command.find("-i") == command.end()) throw runtime_error("Not setting input tag -i");
    if(command["-i"].size() == 0) throw runtime_error("Input file name not specified");
    string inputfilename = command["-i"][0];
    if(command.find("-o") == command.end()) throw runtime_error("Not setting input tag -o");
    if(command["-o"].size() == 0) throw runtime_error("Output file name not specified");
    string outputfilename = command["-o"][0];
    json ofj = parseHMAscii(inputfilename);
    ofstream outputfile(outputfilename.data());
    if(!outputfile.is_open()) throw runtime_error("Output file can't be opened");
    outputfile << ofj << endl;
    outputfile.close();
    return 0;
}