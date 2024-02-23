#ifndef ARGPARSE_HPP
#define ARGPARSE_HPP
#include<iostream>
#include<vector>
#include<map>
#include<string>
using namespace std;

map<string, vector<string>> argparse(int argc, char** argv);
#endif