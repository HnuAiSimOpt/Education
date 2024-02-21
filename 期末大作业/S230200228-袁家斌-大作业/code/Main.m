% 有限元计算-四节点壳
% 
clear; close all;
clc;

%% 定义参数
E = 210000;     % 弹性模量
mu = 0.3;       % 泊松比
P = 10;        % 载荷

Lx = 100;       % 长度 
Ly = 30;        % 宽度
Thickness = 1;  % 厚度
h = -0.5;        % 显示应力的位置，大小为（-0.5*Thickness， 0.5*Thickness）

num_x = 100;
num_y = 30;

scaleFactor = 0.1;   % 后处理缩放因子
%% 生产本构矩阵
D = constitutiveMatrix(E, mu);

%% 生成网格
[~, elemnodes, coords] = rectangularMesh(Lx, Ly, num_x, num_y);
% drawingMesh(elemnodes, coords);

%% 计算总自由度
GDof = 6 * size(coords, 1);

%% 计算刚度矩阵
K = globalStiffness2D(GDof, elemnodes, coords, Thickness, D);

%% 边界条件
% 左边固定
fixedNodes = find(coords(:,2)==0);
bc = zeros(1, 5*length(fixedNodes));
bc(1:5:end) = 6 * fixedNodes - 5;
bc(2:5:end) = 6 * fixedNodes - 4;
bc(3:5:end) = 6 * fixedNodes - 3;
bc(4:5:end) = 6 * fixedNodes - 2;
bc(5:5:end) = 6 * fixedNodes - 1;
% 四周铰接固定
% Ty = find(coords(:,2)==0 | coords(:,2)==Lx);
% Tx = find(coords(:,3)==0 | coords(:,3)==Ly);
% Tz = find(coords(:,3)==0 | coords(:,3)==Ly | coords(:,2)==0 | coords(:,2)==Lx);
% bc = [6*Ty-5; 6*Tx-4; 6*Tz-3];


%% 载荷
% load = distributeLoad(GDof, coords, elemnodes, P );

load = sparse(GDof, 1);
loadNodes = find(coords(:,2)==Lx);
% % 拉伸
load(6*loadNodes-4) = P*Ly / num_y;
load(6*loadNodes(1)-4) = P*Ly / num_y / 2;
load(6*loadNodes(end)-4) = P*Ly / num_y / 2;
% 弯曲
load(6*loadNodes(1)-3) = 50;
load(6*loadNodes(end)-3) = -20;
% load(6*loadNodes-3) = P*Ly / num_y / 1e5;
% load(6*loadNodes(1)-3) = P*Ly / num_y / 2e5;
% load(6*loadNodes(end)-3) = P*Ly / num_y / 2e5;

%% 求解
% 位移
disp = solveDisp(K, bc, load);
deltaDisp = zeros(size(coords,1), 6);
for i = 1:6
    deltaDisp(:,i) = disp(i:6:end);
end

% 应力
stress = solveStress(size(coords, 1), elemnodes, coords, D, disp, h);

%% 后处理
coords(:,4) = 0;
scaleFactor = max([Lx, Ly]) * 0.2 / max(max(abs(deltaDisp(:,1:3)))); %自动缩放
% 位移
subplot(2,1,1);
postDisplacement(elemnodes, coords, deltaDisp(:,1:3), scaleFactor, 'M', 0);
% 应力
subplot(2,1,2);
postStress(elemnodes, coords, scaleFactor, deltaDisp(:,1:3), stress, 'S11', h);