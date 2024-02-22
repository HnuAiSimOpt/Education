% 题目： 悬臂梁左端固定，右端部受集中力 F=100kN 作用,长a=2m,厚度t=5mm,宽b=1m,，材料弹性模量 E= 200GPa、 泊松比v= 0.25
a=2;
b=1;
t=5e-3;
E=200e9;
mu=0.25;
F=1e5; %参数输入
h=0.1; %单元边长
[Node_info,Ele_info]=Meshing(a,b,h); %划分网格
[K,D,BB]=Assembly(Node_info,Ele_info,E,mu,t); %求解单刚并进行总刚与“总应变矩阵”的组装
R=Load(Node_info,F); %节点等效载荷
[KK,RR]=BC(Node_info,K,R); %引入边界条件
u=lsqminnorm(KK,RR); %求解节点位移
[sigma_x,sigma_y,sigma_xy]=Stress(BB,Ele_info,D,u); %求单元应力
figure(1);Plot_u(u,Node_info); %绘制位移分布云图
figure(2);Mise_stress(Node_info,sigma_x_node,sigma_y_node,sigma_xy_node); %绘制 Mises 应力云图