%% ���񻮷֣�������α߳��͵�Ԫ�߳�����������ڵ��š��ڵ�����Ľڵ���Ϣ�������Ԫ��š���Ԫ�����ڵ��ŵĵ�Ԫ��Ϣ��
function [Node_info,Ele_info]=Meshing(a,b,h)
axis equal
hold on
%�ڵ���Ϣ
Num_nodes=0;  %�ڵ���
Node_info=[]; %�ڵ�������Ϣ
for i=0:h:a
 for j=0:h:b
 Num_nodes=Num_nodes+1;
 Node_info=[Node_info;Num_nodes,i,j];
 end
end
%��Ԫ��Ϣ
Num_eles=0;  %��Ԫ���
Ele_info=[]; %��Ԫ��Ϣ
%��ֱ�������ν��е�Ԫ����
for i=h:h:a
 for j=h:h:b
 Num_eles=Num_eles+1;
 node_local=[i-h,j-h;
 i,j-h;
 i-h,j];
 node_list=[];
%���ݵ�Ԫ�ڵ������ҵ���Ӧ�Ľڵ���
 for k=1:3
 row_x=find(abs(Node_info(:,2)-node_local(k,1))<1e-6);
 row_y=find(abs(Node_info(:,3)-node_local(k,2))<1e-6);
 num_node=intersect(row_x,row_y);
 node_list=[node_list;num_node];
 end 
 Ele_info=[Ele_info;Num_eles,node_list'];
%���Ͻ�Ϊֱ�ǵĵ�Ԫ��Ϣ
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

