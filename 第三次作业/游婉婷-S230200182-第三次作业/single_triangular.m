function [Ke,matmtx,kinmtx2]=single_triangular(node_matrix,E,poisson,h)
C=[1,node_matrix(1,2),node_matrix(1,3);
 1,node_matrix(2,2),node_matrix(2,3);
 1,node_matrix(3,2),node_matrix(3,3)];
A=0.5*det(C);  %得到网格单元三角形面积
%得到应变矩阵
kinmtx2=0.5/A*[node_matrix(2,3)-node_matrix(3,3),0,node_matrix(3,3)-node_matrix(1,3),0,node_matrix(1,3)-node_matrix(2,3),0;
 0,node_matrix(3,2)-node_matrix(2,2),0,node_matrix(1,2)-node_matrix(3,2),0,node_matrix(2,2)-node_matrix(1,2);
 node_matrix(3,2)-node_matrix(2,2),node_matrix(2,3)-node_matrix(3,3),node_matrix(1,2)-node_matrix(3,2),...
 node_matrix(3,3)-node_matrix(1,3),node_matrix(2,2)-node_matrix(1,2),node_matrix(1,3)-node_matrix(2,3)];
matmtx=E/(1-poisson^2)*[1,poisson,0;poisson,1,0;0,0,(1-poisson)/2];  %得到材料矩阵
Ke=kinmtx2'*matmtx*kinmtx2*A*h;  %得到单元刚度矩阵
end