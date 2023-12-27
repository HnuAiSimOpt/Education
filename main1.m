%% S230200171 郑小露
%% 求解问题描述如下：长为3m，宽为1.5m的平板结构，其右端面顶点和底点受竖直向下、大小为2000N的力
%% 左端面固定，弹性模量为2.1e7,泊松比为0.25，板厚为0.01m，求解各节点位移和应力，绘制应变云图，结果保存在txt文件中。
clc
clear;
% %-----主程序
%-------------------------------% 初始化 %------------------------------%
% 导入节点信息和杆单元
clc
clear;
node=[1 0 0 0;
      2 0.5 0 0;
      3 1  0 0;
      4 1.5 0 0;
      5 2 0 0;
      6 2.5 0 0;
      7 3 0 0;
      8 0 0.5 0;
      9 0.5 0.5 0;
      10 1 0.5 0;
      11 1.5 0.5 0;
      12 2 0.5 0;
      13 2.5 0.5 0;
      14 3 0.5 0;
      15 0 1 0;
      16 0.5 1 0;
      17 1 1 0;
      18 1.5 1 0;
      19 2 1 0;
      20 2.5 1 0;
      21 3 1 0;
      22 0 1.5 0;
      23 0.5 1.5 0;
      24 1 1.5 0;
      25 1.5 1.5 0;
      26 2 1.5 0;
      27 2.5 1.5 0;
      28 3 1.5 0;];   %节点信息，第一列为节点编号，2~4列分别为x,y,z方向坐标
    ele=[1 1 2 9;
    2 1 9 8;
    3 2 3 10;
    4 2 10 9;
    5 3 4 11;
    6 3 11 10;
    7 4 5 12;
    8 4 12 11;
    9 5 6 13;
    10 5 13 12;
    11 6 7 14;
    12 6 14 13;
    13 8 9 16;
    14 8 16 15;
    15 9 10 17;
    16 9 17 16;
    17 10 11 18;
    18 10 18 17;
    19 11 12 19;
    20 11 19 18;
    21 12 13 20;
    22 12 20 19;
    23 13 14 21;
    24 13 21 20;
    25 15 16 23;
    26 15 23 22;
    27 16 17 24;
    28 16 24 23;
    29 17 18 25;
    30 17 25 24;
    31 18 19 26;
    32 18 26 25;
    33 19 20 27;
    34 19 27 26;
    35 20 21 28;
    36 20 28 27];
num_ele=size(ele, 1);   


 %---物理参数------------------
E = 2.1e7;               % 弹性模量
t = 0.01;               % 单元厚度
v = 0.25;               % 泊松比

n_ele = length(ele(:, 1));   %单元数
%-----------------------------------------------------------------------%

%-------------------------% 组装整体刚度矩阵 %---------------------------%
%组装总体刚度矩阵
dof = length(node(:, 1))*2;       % 自由度数，梁单元每个节点有3个自由度
                               % (横向位移、扭转角位移、弯曲角位移)
f = ones(dof, 1)*1e8;             % 结构整体外载荷矩阵，整体坐标系下
f_loc = zeros(6, 1);              % 单元外载荷矩阵，局部坐标系下
u = ones(dof, 1)*1e6;             % 位移矩阵
K = zeros(dof);                  % 总体刚度矩阵
stress = zeros(n_ele, 1);         % 单元应力矩阵

for i = 1 : n_ele
    k_ele = TriangleElementStiffness(E, v, t, node(ele(i, 2:4), 2:4));
    K = assemTriangle(K, k_ele, ele(i, 2), ele(i, 3), ele(i, 4));
end
%-----------------------------------------------------------------------%

%---------------------------% 定义边界条件 %-----------------------------%
f = zeros(56, 1) ;
f([13, 14, 55, 56]) = [0, -2000, 0, -2000];
u(1)=0; u(2)=0; u(15)=0; u(16)=0; u(29)=0; u(30)=0; u(43)=0; u(44)=0;
%-----------------------------------------------------------------------%

%-------------------------------% 求解 %--------------------------------%
%求解未知自由度
index = [];           % 未知自由度的索引
p = [];               % 未知自由度对应的节点力矩阵
for i = 1:dof
    if u(i) ~= 0
        index = [index, i];
        p = [p; f(i)];
    end
end
u(index) = K(index, index) \ p;    % 高斯消去
f = K * u;

% 单元应力
stress = zeros(num_ele, 3);
x1 = node(:, 2) + u(1:2:56);
y1 = node(:, 3) + u(2:2:56);
%-----------------------------------------------------------------------%

%------------------------------% 可视化 %-------------------------------%
figure;
for i=1 : n_ele
    u1 = [u(2*ele(i, 2)-1);
        u(2*ele(i, 2));
        u(2*ele(i, 3)-1);
        u(2*ele(i, 3));
        u(2*ele(i, 4)-1);
        u(2*ele(i, 4))];
    stress(i, :) = TriangleElementStress(E, v, node(ele(i ,2:4), 2:3), u1, 1)';   % 单元应力计算
    patch(node(ele(i, 2:4), 2), node(ele(i, 2:4), 3), stress(i, 1), 'FaceColor','flat', 'EdgeColor','k');
end
colormap(jet);  % 使用 jet 颜色图
colorbar;  % 显示颜色条

hold on;
figure;
for i=1 : n_ele
    patch(node(ele(i, 2:4),2), node(ele(i, 2:4),3), ...
        'w', 'FaceColor', 'none', 'LineStyle', '-','EdgeColor', 'r');
    hold on;
    patch(x1(ele(i, 2:4)), y1(ele(i, 2:4)), ...
        'w', 'FaceColor', 'none', 'EdgeColor', 'b');
end
% 在箭头的起点和终点坐标
arrow1_start = [3; 2];  % 第一个箭头起点（x，y）
arrow1_end = [3; 1.5];  % 第一个箭头终点（x，y）

arrow2_start = [3; 0];  % 第二个箭头起点（x，y）
arrow2_end = [3; -0.5];  % 第二个箭头终点（x，y）

% 绘制箭头表示外部力
quiver(arrow1_start(1), arrow1_start(2), arrow1_end(1) - arrow1_start(1), arrow1_end(2) - arrow1_start(2), 0, 'r', 'LineWidth', 0.5, 'MaxHeadSize', 1);
quiver(arrow2_start(1), arrow2_start(2), arrow2_end(1) - arrow2_start(1), arrow2_end(2) - arrow2_start(2), 0, 'r', 'LineWidth', 0.5, 'MaxHeadSize', 1);

%-----------------------------------------------------------------------%
 % 保存位移和应力到txt文件
output_file = 'results.txt';

% 打开文件进行写操作
fid = fopen(output_file, 'w');

% 写入节点位移
fprintf(fid, 'Node Displacements:\n');
fprintf(fid, 'Node\tDisplacement_x\tDisplacement_y\n');
for i = 1:length(u)/2
    fprintf(fid, '%d\t%.6f\t%.6f\n', i, u(2*i-1), u(2*i));
end

% 写入单元应力
fprintf(fid, '\nElement Stresses:\n');
fprintf(fid, 'Element\tStress\n');
for i = 1:num_ele
    fprintf(fid, '%d\t%.6f\n', i, stress(i, 1));
end

% 关闭文件
fclose(fid);

disp(['Results have been saved to ' output_file]);

    
    
    