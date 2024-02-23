%% 读取单元与节点数据并且统计单元和节点数量
ID_Element=xlsread('Element-node.xlsx');%读取单元编号及对应节点
ID_Node=xlsread('Node coordinates.xlsx');%读取节点坐标
Num_Element=size(ID_Element,1);%单元数量
Num_Node=size(ID_Node,1);%节点数量
