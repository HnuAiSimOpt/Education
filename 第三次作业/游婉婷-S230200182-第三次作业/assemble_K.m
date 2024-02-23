function [matmtx,K,B]=assemble_K(node,element,E,poisson,h)
nodes_num=size(node,1);
eles_num=size(element,1);
K=zeros(2*nodes_num);
B=[];
for i=1:eles_num
 node_matrix=[element(i,2),node(element(i,2),2),node(element(i,2),3);
 element(i,3),node(element(i,3),2),node(element(i,3),3);
 element(i,4),node(element(i,4),2),node(element(i,4),3)];
 %单元节点矩阵，其第一列为节点编号，第二列为节点横坐标，第三列为节点纵坐标
 [Ke,matmtx,kinmtx2]=single_triangular(node_matrix,E,poisson,h); %求解单元刚度矩阵
 B=[B;kinmtx2]; 
 j=node_matrix(1,1);
 k=node_matrix(2,1);
 m=node_matrix(3,1);
 num=[2*j-1,2*j,2*k-1,2*k,2*m-1,2*m];
 for n1=1:6
 for n2=1:6
 K(num(n1),num(n2))=K(num(n1),num(n2))+Ke(n1,n2);  %刚度矩阵组装
 end
 end
end 