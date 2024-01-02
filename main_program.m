% 第三次作业：线性三角形单元有限元实现   B230200171 李启明
% 条件：模型，长度0.5m，宽度0.25m，厚度0.025m，右端面受均匀分载荷w=3000 KN/m2
% 左端面固定，弹性模量为2.1e4 Pa，泊松比0.3。
% 分析与网格划分：长度方向8个三角形单元，宽度方向2个三角形单元，
% 左边对3个节点进行固定，右边对三个节点施加外载荷
clc
clear;
node=[1 0 0 0;
      2 0.125 0 0;
      3 0.25  0 0;
      4 0.375 0 0;
      5 0.5 0 0;
      6 0 0.125 0;
      7 0.125 0.125 0;
      8 0.25 0.125 0;
      9 0.375 0.125 0;
      10 0.5 0.125 0;
      11 0 0.25 0;
      12 0.125 0.25 0;
      13 0.25 0.25 0;
      14 0.375 0.25 0;
      15 0.5 0.25 0];   %节点信息，第一列为节点编号，2~4列分别为x,y,z方向坐标
ele=[1 1 2 7;
     2 1 7 6;
     3 2 3 8;
     4 2 8 7;
     5 3 4 9;
     6 3 9 8;
     7 4 5 10;
     8 4 10 9;
     9 6 7 12;
     10 6 12 11;
     11 7 8 13;
     12 7 13 12;
     13 8 9 14;
     14 8 14 13;
     15 9 10 15;
     16 9 15 14];        %单元信息，第一列为单元编号，后面各列为单元上的节点号码
 num_ele=size(ele,1);    %单元数
%  for k = 1:num_ele
%     patch(node(ele(k,2:4),2),node(ele(k,2:4),3),'w','FaceColor','none','LineStyle','-','EdgeColor','r');
%  end
%  axis off;
 %---物理参数------------------
E=2.1e7;             %弹性模量
t=0.025;               %单元厚度
miu=0.3;           %泊松比
%-------------------------------
n_ele=length(ele(:,1));   %单元数

%组装总体刚度矩阵
dof=length(node(:,1))*2;             %自由度数，梁单元每个节点有3个自由度，横向位移、扭转角位移、弯曲角位移
f=ones(dof,1)*1e8;             %结构整体外载荷矩阵，整体坐标系下
f_loc=zeros(6,1);              %单元外载荷矩阵，局部坐标系下
u=ones(dof,1)*1e6;             %位移矩阵
K=zeros(dof);                  %总体刚度矩阵
stress=zeros(n_ele,1);         %单元应力矩阵

for i=1:n_ele
    k_ele=TriangleElementStiffness(E,miu,t,node(ele(i,2:4),2:4));
    K=assemTriangle(K,k_ele,ele(i,2),ele(i,3),ele(i,4));
end

%力边界条件
f(9)=4.6875;       %5节点横向力
f(10)=0;          %5节点垂向力
f(19)=9.375;       %10节点横向力
f(20)=0;          %10节点垂向力
f(29)=4.6875;       %15节点横向力
f(30)=0;          %15节点垂向力
f(3)=0;
f(4)=0;
f(13)=0;
f(14)=0;
f(23)=0;
f(24)=0;
f(5)=0;
f(6)=0;
f(15)=0;
f(16)=0;
f(25)=0;
f(26)=0;
f(7)=0;
f(8)=0;
f(17)=0;
f(18)=0;
f(27)=0;
f(28)=0;
%位移边界条件
u(1)=0;
u(2)=0;
u(11)=0;
u(12)=0;
u(21)=0;
u(22)=0;
%求解未知自由度
index=[];      %未知自由度的索引
p=[];          %未知自由度对应的节点力矩阵
for i=1:dof
    if u(i)~=0
        index=[index,i];
        p=[p;f(i)];
    end
end
u(index)=K(index,index)\p;    %高斯消去
f=K*u;

%单元应力
stress=zeros(num_ele,3);
x1=node(:,2)+u(1:2:30);
y1=node(:,3)+u(2:2:30);
figure;
for i=1:n_ele
    u1=[u(2*ele(i,2)-1);u(2*ele(i,2));u(2*ele(i,3)-1);u(2*ele(i,3));u(2*ele(i,4)-1);u(2*ele(i,4))];
    stress(i,:)=TriangleElementStress(E,miu,node(ele(i,2:4),2:3),u1,1)';   %单元应力计算
    patch(node(ele(i,2:4),2),node(ele(i,2:4),3),stress(i,1));
end
hold on;
figure;
for i=1:n_ele
    patch(node(ele(i,2:4),2),node(ele(i,2:4),3),'w','FaceColor','none','LineStyle','-','EdgeColor','b');
    hold on;
    patch(x1(ele(i,2:4)),y1(ele(i,2:4)),'w','FaceColor','none','EdgeColor','r');
end
function k_ele=TriangleElementStiffness(E,miu,t,node_ele)
%TriangleElementStiffness This function returns the element
% stiffness matrix for a Triangle(CST)
% element with modulus of elasticity E,
% Poission's ratio miu, constant thickness t,
% node_ele the node coordinate of element.
% The size of the element stiffness
% matrix is 6 x 6.
%---------node coordinate------
x1=node_ele(1,1);                
y1=node_ele(1,2);
x2=node_ele(2,1);                
y2=node_ele(2,2);
x3=node_ele(3,1);                
y3=node_ele(3,2);
%-------------------------------
A=(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2;   %单元面积
a1=x2*y3-y2*x3;
a2=y1*x3-x1*y3;
a3=x1*y2-y1*x2;
b1=y2-y3;
b2=y3-y1;
b3=y1-y2;
c1=x3-x2;
c2=x1-x3;
c3=x2-x1;
B=1/2/A*[b1 0 b2 0 b3 0;
         0 c1 0 c2 0 c3;
         c1 b1 c2 b2 c3 b3];
D=E/(1-miu^2)*[1 miu 0;
               miu 1 0;
               0 0 (1-miu)/2];     %strss/strain matrix for plane stress
% D=E/(1-miu^2)*[1 miu 0;
%                miu 1 0;
%                0 0 (1-miu)/2];     %strss/strain matrix for plane strain
k_ele=t*A*B'*D*B;    %的单元刚度矩阵
end
function k_t=assemTriangle(k_t,k_ele,node1,node2,node3)
%assemTriangle This function assembles the element stiffness
% matrix k of the plane Triangle element with nodes
% i and j into the global stiffness matrix K.
% This function returns the global stiffness
% matrix K after the element stiffness matrix
% k is assembled.

d(1:2)=2*node1-1:2*node1;
d(3:4)=2*node2-1:2*node2;
d(5:6)=2*node3-1:2*node3;
for ii=1:6
    for jj=1:6
        k_t(d(ii),d(jj))=k_t(d(ii),d(jj))+k_ele(ii,jj);
    end
end
end
function str=TriangleElementStress(E,miu,node_ele,u1,p)
%TriangleElementStress This function returns the element
% stress matrix for a Triangle
% element with modulus of elasticity E,
% Poission's ratio miu,
% node_ele the node coordinate of element，plane stress or plane strain option.

%---------node coordinate------
x1=node_ele(1,1);                
y1=node_ele(1,2);
x2=node_ele(2,1);                
y2=node_ele(2,2);
x3=node_ele(3,1);                
y3=node_ele(3,2);
%-------------------------------
A=(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2;   %单元面积
a1=x2*y3-y2*x3;
a2=y1*x3-x1*y3;
a3=x1*y2-y1*x2;
b1=y2-y3;
b2=y3-y1;
b3=y1-y2;
c1=x3-x2;
c2=x1-x3;
c3=x2-x1;
B=1/2/A*[b1 0 b2 0 b3 0;
         0 c1 0 c2 0 c3;
         c1 b1 c2 b2 c3 b3];
if p==1
    D=E/(1-miu^2)*[1 miu 0;
                   miu 1 0;
                   0 0 (1-miu)/2];           %strss/strain matrix for plane stress
elseif p==2
    D=E/(1+miu)/(1-2*miu)*[1-miu miu 0;
                           miu 1-miu 0;
                       0 0 (1-2*miu)/2];     %strss/strain matrix for plane strain
end
str=D*B*u1;   %单元应力
end
