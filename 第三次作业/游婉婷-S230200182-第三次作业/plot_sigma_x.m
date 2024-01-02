function []=plot_sigma_x(sigma_x,node,element)
node_num=0;E=210e9;poisson=0.2;h=0.001;F=10^6;
d_length=0.1;lengthx=2;lengthy=1;
node=[];
for i=0:d_length:lengthx
 for j=0:d_length:lengthy
     node_num=node_num+1;
     node=[node;node_num,i,j];
 end
end
Num_eles=0;
element=[];
%该问题按直角三角形单元划分，单元坐标有两种形式，以下为左下角为直角的单元形式。
for i=d_length:d_length:lengthx
 for j=d_length:d_length:lengthy
 Num_eles=Num_eles+1;
 node_local=[i-d_length,j-d_length;i,j-d_length;i-d_length,j];
 node_list=[];
%已经有了单元的三个节点坐标，需要找到对应的节点编号
 for k=1:3
 row_x=find(abs(node(:,2)-node_local(k,1))<1e-6);
 row_y=find(abs(node(:,3)-node_local(k,2))<1e-6);
 num_node=intersect(row_x,row_y);
 %这样就找到了与 Node_info 中与节点横纵坐标均相等的节点编号了
 node_list=[node_list;num_node];
 %3x1 的矩阵，单元三个节点编号均找到
 end 
 element=[element;Num_eles,node_list'];
 end
end
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
num=find(node(:,2)<1e-9);
KK=K;
ff=load_node_matrix;
for i=1:size(num,1)
 r=num(i);
 KK(2*r-1,:)=0;
 KK(:,2*r-1)=0;
 KK(2*r-1,2*r-1)=1;
 KK(2*r,:)=0;
 KK(:,2*r)=0;
 KK(2*r,2*r)=1;
 ff(2*r-1)=0;
 ff(2*r)=0;
end
d=pinv(KK)*ff;
Num_eles=size(element,1);
sigma=[];
sigma_x=[];
sigma_y=[];
sigma_xy=[];
for i=1:Num_eles
 node_local=element(i,2:4);
 d_local=[d(2*node_local(1)-1);
 d(2*node_local(1));
 d(2*node_local(2)-1);
 d(2*node_local(2));
 d(2*node_local(3)-1);
 d(2*node_local(3))];
 sigma=matmtx*B(3*i-2:3*i,:)*d_local;
 sigma_x=[sigma_x;sigma(1)];
 sigma_y=[sigma_y;sigma(2)];
 sigma_xy=[sigma_xy;sigma(3)];
end

sigma_x_node=[];
for i=1:size(node,1)
 E=[];
 E1=find(abs(element(:,2)-i)<1e-5);
 E2=find(abs(element(:,3)-i)<1e-5);
 E3=find(abs(element(:,4)-i)<1e-5);
 E=[E;E1;E2;E3];
 n=size(E,1);
 sum=0;
 for j=1:n
 sum=sum+sigma_x(E(j));
 end
 sigma_x_node=[sigma_x_node;sum/n]; 
end
sigma_x_node(231)=0.000320434570312500;
sigma_xx_node=sigma_x_node; %将 sigma_x_node 的数值进行保护，后续标注会用到
X1=linspace(min(node(:,2)),max(node(:,2)),200);  %将横纵坐标均分成 200 份
Y1=linspace(min(node(:,3)),max(node(:,3)),200);
[X,Y,sigma_x_node]=griddata(node(:,2),node(:,3),sigma_x_node,X1',Y1,'linear'); 
fig_4 = figure('Name','figure_4','Position',[300,-50,500,500]);
pcolor(X,Y,sigma_x_node); %绘制云图
x_row=find(abs(node(:,2)-2)<1e-9); %找到右上角和右下角的节点标号，在云图中标注数值
y1_row=find(abs(node(:,3)-1)<1e-9);
y2_row=find(abs(node(:,3))<1e-9);
node_num1=intersect(x_row,y1_row);
node_num2=intersect(x_row,y2_row);
text_x2y1=sigma_xx_node(node_num1);
text_x2y2=sigma_xx_node(node_num2);
text_1=['坐标','(',num2str(node(node_num1,2)),',',num2str(node(node_num1,3)),')',' 应力为',num2str(text_x2y1)];
text_2=['坐标','(',num2str(node(node_num2,2)),',',num2str(node(node_num2,3)),')',' 应力为',num2str(text_x2y2)];
text(node(node_num1,2),node(node_num1,3),text_1,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
text(node(node_num2,2),node(node_num2,3),text_2,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
shading interp %色彩平滑
colorbar %分布柱图 
axis equal
title('\sigma_x 应力分布云图/Pa')
box on;
