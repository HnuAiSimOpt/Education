function K = K_Assemble(K, k, nodes_index,element_type)
%----------------------------------------------%
%function:
%   Assembling the Stiffness Matrix accordinating to the Node index and
%   Element type
%input:
%   K:
%       System Stiffness Matrix
%   k:
%       Element Stiffenss
%   nodes_index:
%       Index of Element
%       for C3N8:[1,2,3,4,5,6,7,8]
%   element_type:
%       Type of Element. ["C3D4"四面体单元,"D3N8"六面体单元]
%----------------------------------------------%

if element_type == "D3N8"
    for row = 1:8
        row_index = nodes_index(row);%刚度矩阵行节点编号
        for col = 1:8
            col_index = nodes_index(col);%刚度矩阵列节点编号
            K((3*row_index-2):row_index*3,(3*col_index-2):col_index*3) = ...
                K((3*row_index-2):row_index*3,(3*col_index-2):col_index*3)...
                + k((3*row-2):3*row,(3*col-2):3*col);
        end
    end
elseif element_type == "C3D4"
    for row = 1:4
        row_index = nodes_index(row);%刚度矩阵行节点编号
        for col = 1:4
            col_index = nodes_index(col);%刚度矩阵列节点编号
            K((3* -2):row_index*3,(3*col_index-2):col_index*3) = ...
                K((3*row_index-2):row_index*3,(3*col_index-2):col_index*3)...
                + k((3*row-2):3*row,(3*col-2):3*col);
        end
    end
end