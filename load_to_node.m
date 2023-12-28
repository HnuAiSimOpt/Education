function load_node_matrix=load_to_node(node,F)
load_node_matrix=zeros(size(node,1)*2,1);
row1_x=find(node(:,2)==2);
row1_y=find(node(:,3)==0);
num1=intersect(row1_x,row1_y);  %找到（2，0）点
row2_x=find(node(:,2)==1);
row2_y=find(node(:,3)==0);
num2=intersect(row2_x,row2_y);  %找到（1，0）点
row3_x=find(node(:,2)==2);
row3_y=find(node(:,3)==1);
num3=intersect(row3_x,row3_y);  %找到（2，1）点
row4_x=find(node(:,2)==1);
row4_y=find(node(:,3)==1);
num4=intersect(row4_x,row4_y);  %找到（1，1）点
load_node_matrix(2*num1)=-F/2;  %给四个节点赋值力
load_node_matrix(2*num3)=-F/2;
load_node_matrix(2*num2)=2*F;
load_node_matrix(2*num4)=-2*F;
end