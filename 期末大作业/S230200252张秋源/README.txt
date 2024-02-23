The program is made by Qiuyuan Zhang, from Hunan University.
这个程序由湖南大学的张秋源制作。

This program is a 4-node tetrahydron FEM solver program, the input and output are json files.
This program is requires several dependencies:

cmake https://cmake.org/ (For building the program)
Eigen https://eigen.tuxfamily.org/
superlu https://github.com/xiaoyeli/superlu
nlohmann/json https://github.com/nlohmann/json

To build this program, you need to
mkdir build
cd build
cmake ..
cmake --build .

There are 2 executable files:
tetra4
The major calculating program, only support tetra4 element, to use it
./tetra4 -i block1-in.json -o block1-out.json

Another program "parsehmscii" is to read hmscii file(And only support tetra4 elements)
./parsehmascii -i block1.hmascii -o block1-in.json

After getting the input json file, you need to mannally set material attribute, like:
"material":{
        "modulus-of-elasticity":2.1e5,
        "poisson-ratio":0.3
    }
in to that input json file.

For costomized input, read the file "element12-input.json" for further details.
The output file contains stress and strain information for each element, as well as the Von
Mises yielding theorm calculated results. And the output file also contains the displace of
each node.

The program contains several major C++ files:
inputjson.cpp (For reading the information of elements and nodes, along with other properties
for this analysis)
tetra4.cpp (Calculating)
tetra4-struct.cpp (Define necessary structs)
resultjson.cpp (Put the calculated result to json file)
