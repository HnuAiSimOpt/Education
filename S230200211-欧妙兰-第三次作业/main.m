%% comment
% purpose：对一矩形薄板，求解其右上角和右下角的位移
% 左端固定，右端部受集中力 F=1000kN 作用；
% 材料弹性模量 E= 210GPa、 泊松比v= 0.2，悬臂梁的厚度(板厚)为 10mm。


clear;
clc;

%% 根据题目要求与计算需要输入参数，主函数计算求解
a=2;
b=1;
t=10e-3;
E=210e9;
mu=0.2;
F=1e6; %参数输入使用“米-帕”单位制
h=0.1; %单元边长（米）,可以控制网格密度
[Node_info,Ele_info]=Meshing(a,b,h); %划分网格

%% 求解单刚并进行总刚与"总应变矩阵"的组装
[K,D,BB]=Assembly(Node_info,Ele_info,E,mu,t); 

%% 节点等效载荷
R=Load(Node_info,F); 
[KK,RR]=BC(Node_info,K,R); %引入边界条件消除总刚奇异性 

%% 求解节点位移
u=lsqminnorm(KK,RR); 
%% 求单元应力
[sigma_x,sigma_y,sigma_xy]=Stress(BB,Ele_info,D,u);

%% 绘制图像结果
figure(1);Plot_u(u,Node_info); %绘制位移分布云图
%figure(2);sigma_x_node=Plot_sx_node(sigma_x,Node_info,Ele_info); %绘制 sigma_x 应力云图，并输出节点平均应力值
%figure(3);sigma_y_node=Plot_sy_node(sigma_y,Node_info,Ele_info); %绘制 sigma_y 应力云图，并输出节点平均应力值
%figure(4);sigma_xy_node=Plot_sxy_node(sigma_xy,Node_info,Ele_info); %绘制 sigma_xy 应力云图，并输出节点平均应力值
%figure(5);Mise_stress(Node_info,sigma_x_node,sigma_y_node,sigma_xy_node); %绘制 Mises 应力云图
