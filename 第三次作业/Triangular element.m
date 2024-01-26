%% 目的：使用线性三角形单元求解求解悬臂薄板的位移，应力。
%%边界条件：矩形薄板右端部受集中力F=1kN作用，材料弹性模量E=10GPa、泊松比v= 0.2，悬臂梁的厚度(板厚)为 10mm
clear all;
clc;

%% 根据题目要求与计算需要输入参数，主函数计算求解
a=1; %矩形板长度(m)
b=0.5; %矩形板宽度(m)
t=10e-3; %板厚(m)
E=10e9; %材料弹性模量(GPa)
mu=0.2; %泊松比
F=1000; %集中力(N)
h=0.1; %单元边长(m)
[Node_info,Ele_info]=Meshing(a,b,h); %划分网格
[K,D,BB]=Assembly(Node_info,Ele_info,E,mu,t); %求解单刚并进行总刚与“总应变矩阵”的组装
R=Load(Node_info,F); %节点等效载荷
[KK,RR]=BC(Node_info,K,R); %引入边界条件消除总刚奇异性 
u=lsqminnorm(KK,RR); %求解节点位移
[sigma_x,sigma_y,sigma_xy]=Stress(BB,Ele_info,D,u); %求解单元应力