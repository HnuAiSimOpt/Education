% 获取网格信息
nodes = mesh.Nodes';
elements = mesh.Elements';

% 定义问题参数
E = 1.0e11;      % 弹性模量 (Pa)
nu = 0.3;       % 泊松比
t = 0.01;        % 材料厚度 (m)


% 计算单元刚度矩阵
num_elements = size(elements, 1);
K_global = zeros(2 * size(nodes,1));

for e = 1:num_elements
    nodes_e = elements(e, :);
    x_e = nodes(nodes_e,1);
    y_e = nodes(nodes_e,2);

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

function [Ke, fe] = linear_triangle_element_stiffness(E, nu, t, x, y)
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