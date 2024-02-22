% 应用边界条件
fixed_nodes = find(nodes(:, 1) == 0 );%左端节点固定
fixed_dofs = [2 * fixed_nodes - 1; 2 * fixed_nodes];
free_dofs = setdiff(1:2 * size(nodes,1), fixed_dofs);

% 施加荷载
force = zeros(2 * size(nodes,1), 1);
force(2 * find(nodes(:, 1) == L)-1, 1) = 1.0e5;  % 在右侧施加均布荷载

% 求解位移场
displacements = zeros(2 * size(nodes,1), 1);
displacements(free_dofs) = K_global(free_dofs, free_dofs) \ force(free_dofs);

nodes_new=nodes;
for i=1:size(nodes,1)
    nodes_new(i,1)=nodes_new(i,1)+displacements(2*i-1);
    nodes_new(i,2)=nodes_new(i,2)+displacements(2*i);
end

% 提取 x 和 y 坐标
x_coords = nodes_new(:, 1);
y_coords = nodes_new(:, 2);

% 绘制散点图
scatter(x_coords, y_coords, 'o', 'filled');

% 添加标签和标题
xlabel('X 坐标');
ylabel('Y 坐标');
title('节点坐标散点图');

% 显示网格
grid on;