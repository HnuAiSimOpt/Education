%% 单元刚度矩阵计算（输入节点、单元信息，输出单元刚度、弹性矩阵，应变矩阵）
function [ke,D,B]=Ke(node_info,E,mu,t)
C=[1,node_info(1,2),node_info(1,3);
 1,node_info(2,2),node_info(2,3);
 1,node_info(3,2),node_info(3,3)];
A=0.5*det(C);
B=0.5/A*[node_info(2,3)-node_info(3,3),0,node_info(3,3)-node_info(1,3),0,node_info(1,3)-node_info(2,3),0;
 0,node_info(3,2)-node_info(2,2),0,node_info(1,2)-node_info(3,2),0,node_info(2,2)-node_info(1,2);
 node_info(3,2)-node_info(2,2),node_info(2,3)-node_info(3,3),node_info(1,2)-node_info(3,2),...
 node_info(3,3)-node_info(1,3),node_info(2,2)-node_info(1,2),node_info(1,3)-node_info(2,3)];
D=E/(1-mu^2)*[1,mu,0;
 mu,1,0;
 0,0,(1-mu)/2];
ke=B'*D*B*A*t;
