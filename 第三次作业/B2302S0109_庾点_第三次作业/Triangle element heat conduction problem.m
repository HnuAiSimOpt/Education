% 有限元线性三角形程序
% 学号：B2302S0109
% 姓名：庾点
%
% 该程序用于求解热传导问题，采用线性三角形元素进行离散。
% 边界条件：左边界温度为0，右边界温度为100，下边界温度为0。
% 施加的边界力：在节点5施加大小为10的力。
%
% 使用的单元类型：线性三角形元素

% 定义网格
x = [0 1 1 0.5 0.5]; % x坐标
y = [0 0 1 0.5 0]; % y坐标

% 定义节点和单元
nodes = [x' y']; % 节点坐标
elements = [1 2 4; 2 3 4; 3 5 4]; % 单元节点索引

% 定义边界条件
boundaryNodes = [1 2 3]; % 边界节点索引
boundaryValues = [0 100 0]; % 边界节点对应的温度值

% 定义材料特性
k = 0.5; % 热导率

% 定义施加力
forceNode = 5; % 施加力的节点索引
forceMagnitude = 10; % 施加力的大小

% 创建稀疏矩阵
numNodes = size(nodes, 1);
A = sparse(numNodes, numNodes);
b = zeros(numNodes, 1);

% 循环遍历每个单元
numElements = size(elements, 1);
for i = 1:numElements
    % 获取当前单元的节点索引
    nodeIndices = elements(i, :);
    
    % 计算当前单元的刚度矩阵
    [kMatrix, fVector] = computeElementStiffness(nodes(nodeIndices, :), k);
    
    % 将单元刚度矩阵和载荷向量添加到总刚度矩阵和载荷向量中
    A(nodeIndices, nodeIndices) = A(nodeIndices, nodeIndices) + kMatrix;
    b(nodeIndices) = b(nodeIndices) + fVector;
end

% 处理边界条件
A(boundaryNodes, :) = 0;
A(sub2ind(size(A), boundaryNodes, boundaryNodes)) = 1;
b(boundaryNodes) = boundaryValues;

% 施加力
b(forceNode) = b(forceNode) + forceMagnitude;

% 解线性方程组
temperatures = A\b;

% 输出结果
disp('节点温度：');
disp(temperatures);

% 计算单元刚度矩阵的函数
function [kMatrix, fVector] = computeElementStiffness(nodes, k)
    % 输入参数：
    % nodes: 单元节点坐标矩阵，每行表示一个节点的坐标 [x1, y1; x2, y2; x3, y3]
    % k: 材料热导率
    
    % 计算单元面积
    A = abs((nodes(2,1)-nodes(1,1))*(nodes(3,2)-nodes(1,2)) - (nodes(3,1)-nodes(1,1))*(nodes(2,2)-nodes(1,2))) / 2;
    
    % 计算单元刚度矩阵
    kMatrix = (k / (4 * A)) * [2, -1, -1; -1, 2, -1; -1, -1, 2];
    
    % 初始化单元载荷向量为零向量
    fVector = zeros(3, 1);
end