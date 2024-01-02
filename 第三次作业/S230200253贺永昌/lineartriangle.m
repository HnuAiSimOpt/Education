% 学号: S230200253
% 姓名: 贺永昌
% 程序的用途: 该程序用于如何使用有限元方法对线性三角形元素进行分析
% 采用了三角形单元，定义了两个带有边界条件的节
% 定义三角形网格的节点坐标（顺序为逆时针）
nodes = [0, 0; 1, 0; 0, 1];

% 定义三角形网格的单元（节点顺序使用 1-based index）
elements = [1, 2, 3];

% 定义边界条件
bc_node = [1, 2]; % 带有边界条件的节点
bc_val = [0, 0]; % 边界条件数值

% 创建刚度矩阵和载荷向量
num_nodes = size(nodes, 1);
stiffness = zeros(num_nodes, num_nodes);
load_vector = zeros(num_nodes, 1);

% 遍历单元计算刚度矩阵和载荷向量
num_elements = size(elements, 1);
for i = 1:num_elements
    % 获取单元节点索引
    node_ids = elements(i, :);
    % 计算单元刚度矩阵
    A = [1, nodes(node_ids(1), :); 1, nodes(node_ids(2), :); 1, nodes(node_ids(3), :)];
    area = 0.5 * det(A);
    B = (1 / (2 * area)) * [nodes(node_ids(2), 2) - nodes(node_ids(3), 2), nodes(node_ids(3), 2) - nodes(node_ids(1), 2), nodes(node_ids(1), 2) - nodes(node_ids(2), 2);
                            nodes(node_ids(3), 1) - nodes(node_ids(2), 1), nodes(node_ids(1), 1) - nodes(node_ids(3), 1), nodes(node_ids(2), 1) - nodes(node_ids(1), 1)];
    stiffness(node_ids, node_ids) = stiffness(node_ids, node_ids) + area * B' * B;
    % 计算载荷向量
    load_vector(node_ids) = load_vector(node_ids) + area / 3;
end

% 处理边界条件
for i = 1:length(bc_node)
    node_id = bc_node(i);
    % 将刚度矩阵的对应行和列置为0
    stiffness(node_id, :) = 0;
    stiffness(node_id, node_id) = 1;
    % 更新载荷向量
    load_vector(node_id) = bc_val(i);
end

% 解线性方程组，得到节点位移
displacements = stiffness \ load_vector;

% 输出结果
disp('节点位移:');
disp(displacements);
