/*Program made by Qiuyuan Zhang*/

#include<iostream>
#include<vector>
#include<map>
#include<string>
using namespace std;

map<string, vector<string>> argparse(int argc, char** argv){
    map<string, vector<string> > m1;
    if(argc > 1){
        vector<string> valstr;
        string minus_tag_str = "begin";
        int i;
        for(i = 1; i < argc; i++){
            string temp = argv[i];
            if(temp[0] == '-') break;
            valstr.push_back(temp);
        }
        //m1["begin"] = valstr;
        //bool minus_tag = false;
        
        while(i < argc){
            string temp = argv[i];
            if(temp[0] == '-'){
                m1[minus_tag_str] = valstr;
                valstr.clear();
                minus_tag_str = temp;
            }
            else{
                valstr.push_back(temp);
            }
            i++;
        }
        m1[minus_tag_str] = valstr;
    }
    return m1;
}


