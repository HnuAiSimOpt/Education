%姓名：李伟
%学号：S230200170


% 定义问题参数
E = 1.0e11;      % 弹性模量 (Pa)
nu = 0.3;       % 泊松比
t = 1.0;        % 材料厚度 (m)

% 定义几何参数
L = 1.0;       % 结构长度 (m)
H = 0.5;        % 结构高度 (m)

% 定义网格参数
num_nodes_x = 20;    % X方向节点数
num_nodes_y = 10;    % Y方向节点数

% 生成节点坐标
x = linspace(0, L, num_nodes_x);
y = linspace(0, H, num_nodes_y);

[X, Y] = meshgrid(x, y);
nodes = [X(:), Y(:)];

% 生成单元连接关系
elements = delaunay(X(:), Y(:));

% 计算单元刚度矩阵
num_elements = size(elements, 1);
K_global = zeros(2 * num_nodes_x * num_nodes_y);

for e = 1:num_elements
    nodes_e = elements(e, :);
    x_e = X(nodes_e);
    y_e = Y(nodes_e);

    [Ke, ~] = linear_triangle_element_stiffness(E, nu, t, x_e, y_e);
    
   % 将单元刚度矩阵添加到全局刚度矩阵
   for i = 1:3
    for j = 1:3
        node_i = nodes_e(i);
        node_j = nodes_e(j);
        K_global(2*node_i-1:2*node_i, 2*node_j-1:2*node_j) = ...
            K_global(2*node_i-1:2*node_i, 2*node_j-1:2*node_j) + Ke(2*i-1:2*i, 2*j-1:2*j);
    end
   end

end

% 应用边界条件
fixed_nodes = find(nodes(:, 1) == 0 );%左端节点固定
fixed_dofs = [2 * fixed_nodes - 1; 2 * fixed_nodes];
free_dofs = setdiff(1:2 * num_nodes_x * num_nodes_y, fixed_dofs);

% 施加荷载
force = zeros(2 * num_nodes_x * num_nodes_y, 1);
force(2 * find(nodes(:, 1) == L), 1) = -1.0e7;  % 在右侧施加均布荷载

% 求解位移场
displacements = zeros(2 * num_nodes_x * num_nodes_y, 1);
displacements(free_dofs) = K_global(free_dofs, free_dofs) \ force(free_dofs);

% 计算位移后的节点坐标
X_displaced = X(:) + displacements(1:2:end);
Y_displaced = Y(:) + displacements(2:2:end);

% 计算坐标轴范围
x_min = min([X(:); X_displaced]);
x_max = max([X(:); X_displaced]);
y_min = min([Y(:); Y_displaced]);
y_max = max([Y(:); Y_displaced]);

% 创建一个新的图形窗口
figure;

% 子图1：位移前的结构
subplot(1, 2, 1);  % 1行2列的子图，这是第一个
plot(X(:), Y(:), 'ko', 'MarkerFaceColor', 'cyan'); % 使用圆圈标记原始位置
hold on;
for i = 1:size(elements, 1)
    % 获取当前三角形的节点索引
    tri_indices = elements(i, :);

    % 获取这些节点的坐标
    tri_x = X(tri_indices);
    tri_y = Y(tri_indices);

    % 添加第一个节点到末尾，以闭合三角形
    tri_x(end+1) = tri_x(1);
    tri_y(end+1) = tri_y(1);

    % 绘制三角形
    plot(tri_x, tri_y, 'k-'); % 使用黑色线条连接节点
    hold on; % 保持图形，以便绘制下一个三角形
end
title('Structure Before Displacement');
xlabel('X');
ylabel('Y');
xlim([x_min, x_max]);
ylim([y_min, y_max]);
hold off;

% 子图2：位移后的结构
subplot(1, 2, 2);  % 1行2列的子图，这是第二个
plot(X_displaced, Y_displaced, 'ko', 'MarkerFaceColor', 'red'); % 使用圆圈标记位移后位置
hold on;
for i = 1:size(elements, 1)
    % 获取当前三角形的节点索引
    tri_indices = elements(i, :);

    % 获取这些节点的位移后坐标
    tri_x_displaced = X_displaced(tri_indices);
    tri_y_displaced = Y_displaced(tri_indices);

    % 添加第一个节点到末尾，以闭合三角形
    tri_x_displaced(end+1) = tri_x_displaced(1);
    tri_y_displaced(end+1) = tri_y_displaced(1);

    % 绘制位移后的三角形
    plot(tri_x_displaced, tri_y_displaced, 'r--'); % 使用红色虚线连接节点
    hold on; % 保持图形，以便绘制下一个三角形
end
title('Structure After Displacement');
xlabel('X');
ylabel('Y');
xlim([x_min, x_max]);
ylim([y_min, y_max]);
hold off;



function [Ke, fe] = linear_triangle_element_stiffness(E, nu, t, x, y)
    % 计算线性三角形单元的刚度矩阵和载荷向量
    % 输入:
    % E - 弹性模量
    % nu - 泊松比
    % t - 材料厚度
    % x, y - 单元节点坐标
    % 输出:
    % Ke - 单元刚度矩阵
    % fe - 单元载荷向量

    % 计算单元面积
    A = 0.5 * abs(det([1, 1, 1; x(:)'; y(:)']));

   % 计算 b_i 和 c_i
    b = [y(2) - y(3), y(3) - y(1), y(1) - y(2)];
    c = [x(3) - x(2), x(1) - x(3), x(2) - x(1)];

    % 计算 B 矩阵
    B = 1 / (2 * A) * [b(1), 0, b(2), 0, b(3), 0;
                       0, c(1), 0, c(2), 0, c(3);
                       c(1), b(1), c(2), b(2), c(3), b(3)];

    % 计算D矩阵
    D = (E / (1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];

    % 计算单元刚度矩阵
    Ke = A * t * B' * D * B;

    % 计算单元载荷向量
    fe = zeros(6, 1);
end
