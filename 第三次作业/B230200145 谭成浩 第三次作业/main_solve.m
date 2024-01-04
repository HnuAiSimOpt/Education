clc
clear;

%线性三角形单元
%谭成浩 B230200145
%一块平板，长为4m，宽为3m，厚度为0.025m，弹性模量为21MPa，泊松比为0.3，
%采用平面线性三角形单元进行应力应变计算。

%节点位置信息
node=[1 0 0 0;
      2 1 0 0;
      3 2  0 0;
      4 3 0 0;
      5 4 0 0;
      6 0 1 0;
      7 1 1 0;
      8 2 1 0;
      9 3 1 0;
      10 4 1 0;
      11 0 2 0;
      12 1 2 0;
      13 2 2 0;
      14 3 2 0;
      15 4 2 0;
      16 0 3 0;
      17 1 3 0;
      18 2 3 0;
      19 3 3 0;
      20 4 3 0];
  
%三角形单元各节点编号
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
     16 9 15 14;
     17 11 12 17;
     18 11 17 16;
     19 12 13 18;
     20 12 18 17;
     21 13 14 19;
     22 13 18 19;
     23 14 15 20;
     24 14 20 19];
 num_ele=size(ele,1);

%材料属性
E=2.1e7;
t=0.025;
miu=0.3; 





n_ele=length(ele(:,1));


dof=length(node(:,1))*2;
f=zeros(dof,1);
u=ones(dof,1);
K=zeros(dof);

%计算单元刚度与组装
for i=1:n_ele
    k_ele=TriangleElementStiffness(E,miu,t,node(ele(i,2:4),2:4));
    K=assemTriangle(K,k_ele,ele(i,2),ele(i,3),ele(i,4));
end

%载荷位置与大小
f(40)=-1e4;

%边界条件
u(1)=0;
u(2)=0;
u(11)=0;
u(12)=0;
u(21)=0;
u(22)=0;
u(31)=0;
u(32)=0;

index=[];
p=[];
for i=1:dof
    if u(i)~=0
        index=[index,i];
        p=[p;f(i)];
    end
end
u(index)=K(index,index)\p;

%单元应力
stress=zeros(num_ele,3);
x1=node(:,2)+u(1:2:end);
y1=node(:,3)+u(2:2:end);
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

function str=TriangleElementStress(E,miu,node_ele,u1,p)

x1=node_ele(1,1);                
y1=node_ele(1,2);
x2=node_ele(2,1);                
y2=node_ele(2,2);
x3=node_ele(3,1);                
y3=node_ele(3,2);

A=(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2;  %单元面积

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
                   0 0 (1-miu)/2];          
elseif p==2
    D=E/(1+miu)/(1-2*miu)*[1-miu miu 0;
                           miu 1-miu 0;
                       0 0 (1-2*miu)/2];    
end
str=D*B*u1;  %单元应力
end

function k_t=assemTriangle(k_t,k_ele,node1,node2,node3)

d(1:2)=2*node1-1:2*node1;
d(3:4)=2*node2-1:2*node2;
d(5:6)=2*node3-1:2*node3;
for ii=1:6
    for jj=1:6
        k_t(d(ii),d(jj))=k_t(d(ii),d(jj))+k_ele(ii,jj);
    end
end
end

function k_ele=TriangleElementStiffness(E,miu,t,node_ele)

x1=node_ele(1,1);                
y1=node_ele(1,2);
x2=node_ele(2,1);                
y2=node_ele(2,2);
x3=node_ele(3,1);                
y3=node_ele(3,2);

A=(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2;   

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
               0 0 (1-miu)/2];     
k_ele=t*A*B'*D*B;    
end