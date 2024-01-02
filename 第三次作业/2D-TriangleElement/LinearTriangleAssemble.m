function y = LinearTriangleAssemble(K, k, i, j, m)
%LinearTriangleAssemble
%                      函数用来组装三角形单元矩阵
%                      K:待组装矩阵
%                      k:刚度矩阵
%                i, j, m:k刚度矩阵对应的
nodes = [i, j, m];
i = 1;
for node_row = nodes
    j = 1;
    for node_colum = nodes
        K(node_row*2-1, node_colum*2-1) = K(node_row*2-1, node_colum*2-1) + k(i, j);
        K(node_row*2-1, node_colum*2) = K(node_row*2-1, node_colum*2) + k(i, j+1);
        K(node_row*2, node_colum*2-1) = K(node_row*2, node_colum*2-1) + k(i+1, j);
        K(node_row*2, node_colum*2) = K(node_row*2, node_colum*2) + k(i+1, j+1);
        j = j+2;
    end
    i = i+2;
end
y = K;
end