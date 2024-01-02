%% 整体刚度矩阵集成函数 （输入节点、单元信息，输出总纲、修改后的应变矩阵）
function [K,D,BB]=Assembly(Node_info,Ele_info,E,mu,t)
Num_nodes=size(Node_info,1);
Num_eles=size(Ele_info,1);
K=zeros(2*Num_nodes);
BB=[];
for i=1:Num_eles
 node_info_local=[Ele_info(i,2),Node_info(Ele_info(i,2),2),Node_info(Ele_info(i,2),3);
 Ele_info(i,3),Node_info(Ele_info(i,3),2),Node_info(Ele_info(i,3),3);
 Ele_info(i,4),Node_info(Ele_info(i,4),2),Node_info(Ele_info(i,4),3)];
 %3x3 的矩阵，第一列为节点编号，二、三列为节点横、纵坐标
 [ke,D,B]=Ke(node_info_local,E,mu,t);
 BB=[BB;B]; 
 j=node_info_local(1,1);
 k=node_info_local(2,1);
 m=node_info_local(3,1);
 num=[2*j-1,2*j,2*k-1,2*k,2*m-1,2*m];
 for n1=1:6
 for n2=1:6
 K(num(n1),num(n2))=K(num(n1),num(n2))+ke(n1,n2);
 end
 end
end
end 