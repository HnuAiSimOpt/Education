% 有限元第三次作业
% 许逸驰 S230200195
% 问题描述：一长度为4m，高1m，厚度0.1m的悬臂梁，其杨氏模量为200GPa，泊松比为0.3。
% 左端固定，受到大小为100N/m^2，方向竖直向下的均布载荷。
% 用平面三角形单元进行有限元分析。

% 清空指令
clear all;
close all;
clc;

% 材料参数
E = 200e9;          % 杨氏模量（Pa）
nu = 0.3;           % 泊松比

% 几何参数
L = 4;              % 长度（m）
H = 1;              % 高度（m）
T = 0.1;            % 厚度（m）

% 载荷参数
q = 100;            % 均布载荷（N/m^2）

% 单元参数
n = 20;             % 单元数量
numNodes = n + 1;   % 节点数量

% 创建节点坐标矩阵
x = linspace(0, L, numNodes)';
y = zeros(numNodes, 1);
z = zeros(numNodes, 1);

% 创建单元连接矩阵
connectivity = [(1:numNodes-1)', (2:numNodes)'];

% 创建全局刚度矩阵和载荷向量
K = zeros(numNodes);
F = zeros(numNodes, 1);

% 计算单元刚度矩阵和载荷向量
for i = 1:n
    node1 = connectivity(i, 1);
    node2 = connectivity(i, 2);
    
    % 单元长度和角度
    dx = x(node2) - x(node1);
    dy = y(node2) - y(node1);
    dz = z(node2) - z(node1);
    L_e = sqrt(dx^2 + dy^2 + dz^2);
    cos_theta = dx / L_e;
    sin_theta = dy / L_e;
    tan_theta = sin_theta / cos_theta;
    
    % 单元刚度矩阵
    ke = (E * T / L_e) * [1, -1; -1, 1];
    
    % 单元载荷向量
    fe = (q * T * L_e / 2) * [1; 1];
    
    % 转换为全局坐标系
    Te = [cos_theta, sin_theta; -sin_theta, cos_theta];
    Ke = Te' * ke * Te;
    Fe = Te' * fe;
    
    % 组装全局刚度矩阵和载荷向量
    K(node1:node2, node1:node2) = K(node1:node2, node1:node2) + Ke;
    F(node1:node2) = F(node1:node2) + Fe;
end

% 边界条件（左端固定）
fixedNode = 1;
K(fixedNode, :) = 0;
K(fixedNode, fixedNode) = 1;
F(fixedNode) = 0;

% 解方程，求解位移向量
u = K \ F;

% 计算应力和应变
strain = zeros(n, 1);
stress = zeros(n, 1);
for i = 1:n
    node1 = connectivity(i, 1);
    node2 = connectivity(i, 2);
    
    % 单元长度和角度
    dx = x(node2) - x(node1);
    dy = y(node2) - y(node1);
    dz = z(node2) - z(node1);
    L_e = sqrt(dx^2 + dy^2 + dz^2);
    cos_theta = dx / L_e;
    sin_theta = dy / L_e;
    
    % 单元位移向量
    ue = [u(node1); u(node2)];
    
    % 单元应变
    epsilon_e = (1 / L_e) * [-cos_theta, cos_theta] * ue;
    
    % 单元应力
    sigma_e = E * epsilon_e;
    
    % 存储应变和应力
    strain(i) = epsilon_e;
    stress(i) = sigma_e;
end

% 绘制位移图形
figure;
plot(x, u, 'b-o');
xlabel('x（m）');
ylabel('位移（m）');
title('悬臂梁的位移分布');

% 绘制应变图形
figure;
plot(x(1:end-1), strain, 'r-o');
xlabel('x（m）');
ylabel('应变');
title('悬臂梁的应变分布');

% 绘制应力图形
figure;
plot(x(1:end-1), stress, 'g-o');
xlabel('x（m）');
ylabel('应力（Pa）');
title('悬臂梁的应力分布');
