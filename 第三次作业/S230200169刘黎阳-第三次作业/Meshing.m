%% 网格划分（输入矩形边长和单元边长，输出包含节点编号、节点坐标的节点信息与包含单元编号、单元包含节点编号的单元信息）
function [Node_info,Ele_info]=Meshing(a,b,h)
axis equal
hold on
%节点信息
Num_nodes=0;  %节点编号
Node_info=[]; %节点坐标信息
for i=0:h:a
 for j=0:h:b
 Num_nodes=Num_nodes+1;
 Node_info=[Node_info;Num_nodes,i,j];
 end
end
%单元信息
Num_eles=0;  %单元编号
Ele_info=[]; %单元信息
%以直角三角形进行单元划分
for i=h:h:a
 for j=h:h:b
 Num_eles=Num_eles+1;
 node_local=[i-h,j-h;
 i,j-h;
 i-h,j];
 node_list=[];
%根据单元节点坐标找到对应的节点编号
 for k=1:3
 row_x=find(abs(Node_info(:,2)-node_local(k,1))<1e-6);
 row_y=find(abs(Node_info(:,3)-node_local(k,2))<1e-6);
 num_node=intersect(row_x,row_y);
 node_list=[node_list;num_node];
 end 
 Ele_info=[Ele_info;Num_eles,node_list'];
%右上角为直角的单元信息
 Num_eles=Num_eles+1;
 node_local=[i,j;i-h,j;i,j-h];
 node_list=[];
 for k=1:3
 row_x=find(abs(Node_info(:,2)-node_local(k,1))<1e-6);
 row_y=find(abs(Node_info(:,3)-node_local(k,2))<1e-6);
 num_node=intersect(row_x,row_y);
 node_list=[node_list;num_node];
 end
 Ele_info=[Ele_info;Num_eles,node_list'];
 end

