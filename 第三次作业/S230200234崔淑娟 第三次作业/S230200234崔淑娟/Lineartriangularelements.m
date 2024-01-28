%S230200234 崔淑娟
%问题描述：一矩形薄板长为2m宽为1m,宽边左端固定右端部受集中力 F=1000kN 作用，材料弹性模量 E= 210GPa、 泊松比v = 0.2，悬臂梁的厚度(板厚)为10mm，计算弹性版右端位移与应力并输出应力云图。
clear all;
clc;
 
%% 根据题目要求与计算需要输入参数，主函数计算求解
lengthx=2;%x方向长度
lengthy=1;%y方向长度
t=10e-3; %厚度
emodule=210e9;%弹性模量
poisson=0.2;%泊松比
fload=1e6; %施加力
h=0.1; %单元边长,控制网格密度
[Node_info,Ele_info]=Meshing(lengthx,lengthy,h); %划分网格
[K,D,BB]=Assembly(Node_info,Ele_info,emodule,poisson,t); %求解单刚并进行总刚与“总应变矩阵”的组装
R=Load(Node_info,fload); %节点等效载荷
[KK,RR]=BC(Node_info,K,R); %引入边界条件消除总刚奇异性 
u=lsqminnorm(KK,RR); %求解节点位移
[sigma_x,sigma_y,sigma_xy]=Stress(BB,Ele_info,D,u); %求单元应力
figure(1);Plot_u(u,Node_info); %绘制位移分布云图
figure(2);sigma_x_node=Plot_sx_node(sigma_x,Node_info,Ele_info); %绘制 sigma_x 应力云图，并输出节点平均应力值
figure(3);sigma_y_node=Plot_sy_node(sigma_y,Node_info,Ele_info); %绘制 sigma_y 应力云图，并输出节点平均应力值
figure(4);sigma_xy_node=Plot_sxy_node(sigma_xy,Node_info,Ele_info); %绘制 sigma_xy 应力云图，并输出节点平均应力值
figure(5);Mise_stress(Node_info,sigma_x_node,sigma_y_node,sigma_xy_node); %绘制 Mises 应力云图
 
%% 网格划分函数（输入物体边长和单元边长，将输出包含节点编号、节点坐标的节点信息与包含单元编号、单元包含节点编号的单元信息）
function [Node_info,Ele_info]=Meshing(a,b,h)
axis equal
hold on
%节点信息
Num_nodes=0;
Node_info=[];
for i=0:h:a
    for j=0:h:b
        Num_nodes=Num_nodes+1;
        Node_info=[Node_info;Num_nodes,i,j];
    end
end
%单元信息
Num_eles=0;
Ele_info=[];
%该问题按为左下角为直角的单元形式。
for i=h:h:a
    for j=h:h:b
        Num_eles=Num_eles+1;
        node_local=[i-h,j-h;
            i,j-h;
            i-h,j];
        node_list=[];%找到对应的节点编号
        for k=1:3
            row_x=find(abs(Node_info(:,2)-node_local(k,1))<1e-6);
            row_y=find(abs(Node_info(:,3)-node_local(k,2))<1e-6);
            num_node=intersect(row_x,row_y);
            node_list=[node_list;num_node];
        end
        Ele_info=[Ele_info;Num_eles,node_list'];%右上角为直角的单元信息
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
end
end
 
%% 整体刚度矩阵集成函数 （输入节点、单元信息，输出总纲、修改后的应变矩阵）
function [K,D,BB]=Assembly(Node_info,Ele_info,E,mu,t)
Num_nodes=size(Node_info,1);
Num_eles=size(Ele_info,1);
K=zeros(2*Num_nodes);
BB=[];
for i=1:Num_eles
    node_info_local=[Ele_info(i,2),Node_info(Ele_info(i,2),2),Node_info(Ele_info(i,2),3);
    Ele_info(i,3),Node_info(Ele_info(i,3),2),Node_info(Ele_info(i,3),3);
    Ele_info(i,4),Node_info(Ele_info(i,4),2),Node_info(Ele_info(i,4),3)];%3x3 的矩阵，第一列为节点编号，二、三列为节点横、纵坐标
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
 
%% 单元刚度矩阵计算函数 （输入节点、单元信息，输出单刚、弹性矩阵、应变矩阵）
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
end
 
%% 节点载荷列阵
function R=Load(Node_info,F)
R=zeros(size(Node_info,1)*2,1);
row1_x=find(Node_info(:,2)==2);
row1_y=find(Node_info(:,3)==0);
num1=intersect(row1_x,row1_y);
row2_x=find(Node_info(:,2)==2);
row2_y=find(Node_info(:,3)==1);
num2=intersect(row2_x,row2_y);
R(2*num1)=-F/2;
R(2*num2)=-F/2;
end
 
%% 引入边界条件，修改总刚和载荷列阵 （划一置零法）
function [KK,RR]=BC(Node_info,K,R)
num=find(Node_info(:,2)<1e-9);
KK=K;
RR=R;
for i=1:size(num,1)
    r=num(i);
    KK(2*r-1,:)=0;
    KK(:,2*r-1)=0;
    KK(2*r-1,2*r-1)=1;
    KK(2*r,:)=0;
    KK(:,2*r)=0;
    KK(2*r,2*r)=1;
    RR(2*r-1)=0;
    RR(2*r)=0;
end
end
 
%% 单元应力计算函数 （输入应变总矩阵、单元总信息、计算得到的位移解，输出每个单元的应力）
function [sigma_x,sigma_y,sigma_xy]=Stress(BB,Ele_info,D,u)
Num_eles=size(Ele_info,1);
sigma=[];
sigma_x=[];
sigma_y=[];
sigma_xy=[];
for i=1:Num_eles
    node_local=Ele_info(i,2:4);
    u_local=[u(2*node_local(1)-1);
    u(2*node_local(1));
    u(2*node_local(2)-1);
    u(2*node_local(2));
    u(2*node_local(3)-1);
    u(2*node_local(3))];
    sigma=D*BB(3*i-2:3*i,:)*u_local;
    sigma_x=[sigma_x;sigma(1)];
    sigma_y=[sigma_y;sigma(2)];
    sigma_xy=[sigma_xy;sigma(3)];
end
end
 
%% 后处理函数 绘制节点位移云图
function []=Plot_u(u,Node_info)
u1=[];
u2=[];
for i=1:231
    u1=[u1;u(2*i-1)];
    u2=[u2;u(2*i)];
end
u11=u1;
u22=u2; 
X1=linspace(min(Node_info(:,2)),max(Node_info(:,2)),100); %将坐标均分成 100 份
Y1=linspace(min(Node_info(:,3)),max(Node_info(:,3)),100);
[X,Y,u1]=griddata(Node_info(:,2),Node_info(:,3),u1,X1',Y1,'v4'); %'v4'代表插值方法为 matlab4 样条函数内插
[X,Y,u2]=griddata(Node_info(:,2),Node_info(:,3),u2,X1',Y1,'v4');
subplot(2,1,1);
pcolor(X,Y,u1);
shading interp %色彩平滑
colorbar %分布柱图 
axis equal
title('水平方向位移分布云图')
box on;
x_row=find(abs(Node_info(:,2)-2)<1e-9); %找到右上角和右下角的节点标号，在云图中标注数值
y1_row=find(abs(Node_info(:,3)-1)<1e-9);
y2_row=find(abs(Node_info(:,3))<1e-9);
node_num1=intersect(x_row,y1_row);
node_num2=intersect(x_row,y2_row);
text_x2y1=u11(node_num1);
text_x2y2=u11(node_num2);
text_1=['[','(',num2str(Node_info(node_num1,2)),',',num2str(Node_info(node_num1,3)),')',',',num2str(text_x2y1),']'];
text_2=['[','(',num2str(Node_info(node_num2,2)),',',num2str(Node_info(node_num2,3)),')',',',num2str(text_x2y2),']'];
text(Node_info(node_num1,2),Node_info(node_num1,3),text_1,'HorizontalAlignment', 'right');
text(Node_info(node_num2,2),Node_info(node_num2,3),text_2,'HorizontalAlignment', 'right');
subplot(2,1,2);
pcolor(X,Y,u2);
shading interp %色彩平滑
colorbar %分布柱图 
axis equal
title('竖直方向位移分布云图/m')
box on;
x_row=find(abs(Node_info(:,2)-2)<1e-9); %找到右上角和右下角的节点标号，在云图中标注数值
y1_row=find(abs(Node_info(:,3)-1)<1e-9);
y2_row=find(abs(Node_info(:,3))<1e-9);
node_num1=intersect(x_row,y1_row);
node_num2=intersect(x_row,y2_row);
text_x2y1=u22(node_num1);
text_x2y2=u22(node_num2);
text_1=['[','(',num2str(Node_info(node_num1,2)),',',num2str(Node_info(node_num1,3)),')',',',num2str(text_x2y1),']'];
text_2=['[','(',num2str(Node_info(node_num2,2)),',',num2str(Node_info(node_num2,3)),')',',',num2str(text_x2y2),']'];
text(Node_info(node_num1,2),Node_info(node_num1,3),text_1,'HorizontalAlignment', 'right');
text(Node_info(node_num2,2),Node_info(node_num2,3),text_2,'HorizontalAlignment', 'right');
end
 
%% 后处理函数 绘制 sigma_x 云图
%每一个节点的应力值取周围单元的应力平均值
function sigma_xx_node=Plot_sx_node(sigma_x,Node_info,Ele_info)
sigma_x_node=[];
for i=1:size(Node_info,1)
    E=[];
    E1=find(abs(Ele_info(:,2)-i)<1e-9);
    E2=find(abs(Ele_info(:,3)-i)<1e-9);
    E3=find(abs(Ele_info(:,4)-i)<1e-9);
    E=[E;E1;E2;E3];
    n=size(E,1);
    sx_sum=0;
    for j=1:n
        sx_sum=sx_sum+sigma_x(E(j));
    end
    sigma_x_node=[sigma_x_node;sx_sum/n];
end
sigma_xx_node=sigma_x_node; 
X1=linspace(min(Node_info(:,2)),max(Node_info(:,2)),100); %将坐标均分成 100 份
Y1=linspace(min(Node_info(:,3)),max(Node_info(:,3)),100);
[X,Y,sigma_x_node]=griddata(Node_info(:,2),Node_info(:,3),sigma_x_node,X1',Y1,'v4'); 
pcolor(X,Y,sigma_x_node); %绘制云图
x_row=find(abs(Node_info(:,2)-2)<1e-9); %找到右上角和右下角的节点标号，在云图中标注数值
y1_row=find(abs(Node_info(:,3)-1)<1e-9);
y2_row=find(abs(Node_info(:,3))<1e-9);
node_num1=intersect(x_row,y1_row);
node_num2=intersect(x_row,y2_row);
text_x2y1=sigma_xx_node(node_num1);
text_x2y2=sigma_xx_node(node_num2);
text_1=['[','(',num2str(Node_info(node_num1,2)),',',num2str(Node_info(node_num1,3)),')',',',num2str(text_x2y1),']'];
text_2=['[','(',num2str(Node_info(node_num2,2)),',',num2str(Node_info(node_num2,3)),')',',',num2str(text_x2y2),']'];
text(Node_info(node_num1,2),Node_info(node_num1,3),text_1,'HorizontalAlignment', 'right');
text(Node_info(node_num2,2),Node_info(node_num2,3),text_2,'HorizontalAlignment', 'right');
shading interp %色彩平滑
colorbar %分布柱图 
axis equal
title('\sigma_x 应力分布云图/Pa')
box on;
end
 
%% 后处理函数 绘制 sigma_y 应力云图
function sigma_yy_node=Plot_sy_node(sigma_y,Node_info,Ele_info)
sigma_y_node=[];
for i=1:size(Node_info,1)
    E=[];
    E1=find(abs(Ele_info(:,2)-i)<1e-9);
    E2=find(abs(Ele_info(:,3)-i)<1e-9);
    E3=find(abs(Ele_info(:,4)-i)<1e-9);
    E=[E;E1;E2;E3];
    n=size(E,1);
    sy_sum=0;
    for j=1:n
        sy_sum=sy_sum+sigma_y(E(j));
    end
    sigma_y_node=[sigma_y_node;sy_sum/n];
end
sigma_yy_node=sigma_y_node; %将 sigma_y_node 的数值进行保护，后续标注会用到
X1=linspace(min(Node_info(:,2)),max(Node_info(:,2)),100); %将坐标均分成 100 份
Y1=linspace(min(Node_info(:,3)),max(Node_info(:,3)),100);
[X,Y,sigma_y_node]=griddata(Node_info(:,2),Node_info(:,3),sigma_y_node,X1',Y1,'v4'); 
pcolor(X,Y,sigma_y_node);
x_row=find(abs(Node_info(:,2)-2)<1e-9); %找到右上角和右下角的节点标号，在云图中标注数值
y1_row=find(abs(Node_info(:,3)-1)<1e-9);
y2_row=find(abs(Node_info(:,3))<1e-9);
node_num1=intersect(x_row,y1_row);
node_num2=intersect(x_row,y2_row);
text_x2y1=sigma_yy_node(node_num1);
text_x2y2=sigma_yy_node(node_num2);
text_1=['[','(',num2str(Node_info(node_num1,2)),',',num2str(Node_info(node_num1,3)),')',',',num2str(text_x2y1),']'];
text_2=['[','(',num2str(Node_info(node_num2,2)),',',num2str(Node_info(node_num2,3)),')',',',num2str(text_x2y2),']'];
text(Node_info(node_num1,2),Node_info(node_num1,3),text_1,'HorizontalAlignment', 'right');
text(Node_info(node_num2,2),Node_info(node_num2,3),text_2,'HorizontalAlignment', 'right');
shading interp %色彩平滑
colorbar %分布柱图 
axis equal
title('\sigma_y 应力分布云图/Pa')
box on;
end
 
%% 后处理函数 绘制剪应力云图
function sigma_xxyy_node=Plot_sxy_node(sigma_xy,Node_info,Ele_info)
sigma_xy_node=[];
for i=1:size(Node_info,1)
    E=[];
    E1=find(abs(Ele_info(:,2)-i)<1e-9);
    E2=find(abs(Ele_info(:,3)-i)<1e-9);
    E3=find(abs(Ele_info(:,4)-i)<1e-9);
    E=[E;E1;E2;E3];
    n=size(E,1);
    sxy_sum=0;
    for j=1:n
        sxy_sum=sxy_sum+sigma_xy(E(j));
    end
    sigma_xy_node=[sigma_xy_node;sxy_sum/n];
end
sigma_xxyy_node=sigma_xy_node; %将 sigma_xy_node 的数值进行保护，后续标注会用到
X1=linspace(min(Node_info(:,2)),max(Node_info(:,2)),100); %将坐标均分成 100 份
Y1=linspace(min(Node_info(:,3)),max(Node_info(:,3)),100);
[X,Y,sigma_xy_node]=griddata(Node_info(:,2),Node_info(:,3),sigma_xy_node,X1',Y1,'v4'); 
pcolor(X,Y,sigma_xy_node);
x_row=find(abs(Node_info(:,2)-2)<1e-9); %找到右上角和右下角的节点标号，在云图中标注数值
y1_row=find(abs(Node_info(:,3)-1)<1e-9);
y2_row=find(abs(Node_info(:,3))<1e-9);
node_num1=intersect(x_row,y1_row);
node_num2=intersect(x_row,y2_row);
text_x2y1=sigma_xxyy_node(node_num1);
text_x2y2=sigma_xxyy_node(node_num2);
text_1=['[','(',num2str(Node_info(node_num1,2)),',',num2str(Node_info(node_num1,3)),')',',',num2str(text_x2y1),']'];
text_2=['[','(',num2str(Node_info(node_num2,2)),',',num2str(Node_info(node_num2,3)),')',',',num2str(text_x2y2),']'];
text(Node_info(node_num1,2),Node_info(node_num1,3),text_1,'HorizontalAlignment', 'right');
text(Node_info(node_num2,2),Node_info(node_num2,3),text_2,'HorizontalAlignment', 'right');
shading interp %色彩平滑
colorbar %分布柱图 
axis equal
title('\tau_{xy}分布云图/Pa')
box on;
end
 
%% 后处理函数 绘制 Mises 应力分布云图
function []=Mise_stress(Node_info,sigma_xx_node,sigma_yy_node,sigma_xxyy_node)% 先求主应力
sigma_m=[];
for i=1:size(Node_info,1)
    sigma1=(sigma_xx_node(i)+sigma_yy_node(i))/2+0.5*sqrt((sigma_xx_node(i)-sigma_yy_node(i))^2+4*sigma_xxyy_node(i)*sigma_xxyy_node(i));
    sigma2=(sigma_xx_node(i)+sigma_yy_node(i))/2-0.5*sqrt((sigma_xx_node(i)-sigma_yy_node(i))^2+4*sigma_xxyy_node(i)*sigma_xxyy_node(i));% 求 Von-Mise 等效应力
    sigma_m=[sigma_m;sqrt(sigma1^2+sigma2^2-sigma1*sigma2)];
end
sigma_mm=sigma_m; %将 sigma_m 的数值进行保护，后续标注会用到
X1=linspace(min(Node_info(:,2)),max(Node_info(:,2)),100); %将坐标均分成 100 份
Y1=linspace(min(Node_info(:,3)),max(Node_info(:,3)),100);
[X,Y,sigma_m]=griddata(Node_info(:,2),Node_info(:,3),sigma_m,X1',Y1,'v4'); 
pcolor(X,Y,sigma_m); %绘制云图
x_row=find(abs(Node_info(:,2)-2)<1e-9); %根据坐标找到右上角和右下角的节点标号，在云图中标注数值
y1_row=find(abs(Node_info(:,3)-1)<1e-9);
y2_row=find(abs(Node_info(:,3))<1e-9);
node_num1=intersect(x_row,y1_row);
node_num2=intersect(x_row,y2_row);
text_x2y1=sigma_mm(node_num1);
text_x2y2=sigma_mm(node_num2);
text_1=['[','(',num2str(Node_info(node_num1,2)),',',num2str(Node_info(node_num1,3)),')',',',num2str(text_x2y1),']'];
text_2=['[','(',num2str(Node_info(node_num2,2)),',',num2str(Node_info(node_num2,3)),')',',',num2str(text_x2y2),']'];
text(Node_info(node_num1,2),Node_info(node_num1,3),text_1,'HorizontalAlignment', 'right');
text(Node_info(node_num2,2),Node_info(node_num2,3),text_2,'HorizontalAlignment', 'right');
shading interp %色彩平滑
colorbar %分布柱图 
axis equal
title('Mises 应力分布云图/Pa')
box on;
end