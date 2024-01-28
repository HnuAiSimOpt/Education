% 张晓童
% S230200222
% 对一个长4m，高0.8m，厚0.2m的悬臂梁进行分析，并绘制结构的形变和应力分布图。
% 该悬臂梁左侧固定，受到大小为100N/㎡，方向为y轴负方向的面力。
% 在该程序中，采用了平面直角三角形单元进行有限元分析。

clear;%清除工作区
clc;%清楚命令行
close all;%关闭所有Figure窗口

E=2e11;%弹性模量，Pa
u=0.4;%泊松比
q=100;%外载荷（面力），N/㎡
t=0.2;%厚度，m
D=[E/(1-u^2),u*E/(1-u^2),0;u*E/(1-u^2),E/(1-u^2),0;0,0,E/(2*(1+u))];%应力应变矩阵
bt=(3+0.5407762)*10^11;%乘大数法的大数
long=4;%长6m
height=0.8;%高1m
ele_size=[0.1,0.1];%单元的两个直角边尺寸
ele_H=long/ele_size(1);%水平单元数，40个
ele_V=height/ele_size(2);%竖向单元数，8个
ele_num=ele_H*ele_V*2;%单元数，640个
node_num=(ele_H+1)*(ele_V+1);%节点数，369个
DOF=2;%每个节点自由度
A=1/2*ele_size(1)*ele_size(2);%单元面积，0.005

%定义每个单元的节点号
num=0;
for i=1:ele_H
    for j=1:ele_V
        num=num+1;
        ele(num,1)=1+(j-1)+(i-1)*(ele_V+1); 
        ele(num,2)=1+(j-1)+i*(ele_V+1);
        ele(num,3)=2+(j-1)+(i-1)*(ele_V+1);
        num=num+1;
        ele(num,1)=1+j+(i-1)*(ele_V+1);
        ele(num,2)=1+(j-1)+i*(ele_V+1);
        ele(num,3)=2+(j-1)+i*(ele_V+1);
    end
end

%定义节点坐标
num=0;
for i=1:ele_H+1
    for j=1:ele_V+1
        num=num+1;
        node(num,1)=(i-1)*ele_size(1);%x
        node(num,2)=(j-1)*ele_size(2);%y
    end
end
        
%每个单元自由度在整体矩阵中的索引位置
ele_DOF=zeros(ele_num,DOF*3);
for i=1:ele_num
    for j=1:length(ele(i,:))%三角形三个节点循环
        ele_DOF(i,1+2*(j-1):2+2*(j-1))=2*(ele(i,j)-1)+1:2*(ele(i,j)-1)+2;
    end
end

%组装系统刚度矩阵
KZ=zeros(node_num*2,node_num*2);
for i=1:ele_num 
    ke=ele_stiff_matrix(ele,node,i,A,D,t);%单元刚度矩阵
    KZ(ele_DOF(i,:),ele_DOF(i,:)) = KZ(ele_DOF(i,:),ele_DOF(i,:))+ke;
end

%计算总节点载荷阵列
P=[]; 
P(2*(ele_V+1))=-1/2*q*t*ele_size(1);
for i=2:ele_H
      P(18+2*(ele_V+1)*(i-1))=-q*t*ele_size(1);
end
P(node_num*2)=-1/2*q*t*ele_size(1);
PZ=P';%总节点载荷列阵

%乘大数法施加约束条件，最左侧节点全部固定
cons_DOF=(1:(ele_V+1)*2);%约束自由度编号
Disp=zeros(length(cons_DOF),1);%位移为0
for i=1:length(cons_DOF)
      KZ(cons_DOF(i),cons_DOF(i))=bt*KZ(cons_DOF(i),cons_DOF(i));
      PZ(cons_DOF(i))=Disp(i)*KZ(cons_DOF(i),cons_DOF(i));
end

a=inv(KZ)*PZ;%结点位移列阵
for i=1:ele_num %640个单元
      sgm(:,i)=stress_calculate(ele,node,i,A,D,t,a);%单元应力
end
  
%重新排列位移矩阵
SF=1/10*1/max(abs(a));%位移放大系数
new_node=node+SF*reshape(a,2,length(a)/2)';%计算变形后的坐标

%形变
figure
patch('Faces',ele,'Vertices',node,'FaceColor','none','EdgeColor','g')
hold on;
patch('Faces',ele,'Vertices',new_node,'FaceColor','none')
xlabel('x')
ylabel('y')
title('形变')
axis equal

%应力  
figure
patch('Faces',ele,'Vertices',node,'FaceVertexCData',sgm(1,:)','FaceColor','flat')
axis equal
xlabel('x')
ylabel('y')
title('应力\sigma  Pa')
colorbar
colormap(jet)
