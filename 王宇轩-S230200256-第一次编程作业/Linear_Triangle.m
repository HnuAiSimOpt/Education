% 学号: S230200256
% 姓名: 王宇轩
% 程序的用途: 该程序用于如何使用有限元方法对线性三角形元素进行分析:
% 定义了一个四个节点的三角形板，施加了在节点4上的负y方向力，以及节点1, 2和3的固定边界条件。


nodes = [0, 0; 1, 0; 0, 1; 1, 1]; % 节点坐标
elements = [1 2 3; 2 3 4]; % 单元节点索引
loads = [4, 0, -1]; % 施加的力
constraints = [1, 0, 0; 2, 0, 0; 3, 0, 0]; % 边界条件
E = 210e9; % 弹性模量
nu = 0.3; % 泊松比
thickness = 0.01; % 板的厚度

u = fem_analysis(nodes, elements, loads, constraints, E, nu, thickness);
disp(u);


function [u] = fem_analysis(nodes, elements, loads, constraints, E, nu, thickness)
    % nodes: 节点坐标数组 [x1, y1; x2, y2; ...]
    % elements: 单元节点索引数组 [n1, n2, n3; ...]，表示第i个单元连接了节点n1, n2和n3
    % loads: 节点施加的力数组 [n, fx, fy; ...]，表示施加在节点n上的力的x和y分量
    % constraints: 节点边界条件数组 [n, u1, u2; ...]，表示节点n的x和y方向的位移
    % E: 弹性模量
    % nu: 泊松比
    % thickness: 板的厚度
    
    numNodes = size(nodes, 1); % 节点数
    numElements = size(elements, 1); % 单元数
    
    % 计算单元刚度矩阵和全局刚度矩阵
    globalK = zeros(2 * numNodes); % 全局刚度矩阵
    for i = 1:numElements
        elementNodes = elements(i, :);
        elementCoordinates = nodes(elementNodes, :);
        elementK = compute_element_stiffness(elementCoordinates, E, nu, thickness);
        globalK = assemble_element_stiffness(globalK, elementK, elementNodes);
    end
    
    % 处理施加的力和节点边界条件
    globalF = zeros(2 * numNodes, 1); % 全局力向量
    for i = 1:size(loads, 1)
        node = loads(i, 1);
        fx = loads(i, 2);
        fy = loads(i, 3);
        globalF(2 * node - 1) = fx;
        globalF(2 * node) = fy;
    end
    
    globalU = globalK \ globalF; % 计算位移向量
    
    % 利用位移向量计算节点位移数组
    u = reshape(globalU, 2, [])';
    
    % 处理节点边界条件
    for i = 1:size(constraints, 1)
        node = constraints(i, 1);
        u(node, :) = [constraints(i, 2), constraints(i, 3)];
    end
end

function [elementK] = compute_element_stiffness(elementCoordinates, E, nu, thickness)
    % 计算单元刚度矩阵
    x1 = elementCoordinates(1, 1);
    y1 = elementCoordinates(1, 2);
    x2 = elementCoordinates(2, 1);
    y2 = elementCoordinates(2, 2);
    x3 = elementCoordinates(3, 1);
    y3 = elementCoordinates(3, 2);
    
    area = abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / 2;
    D = E * thickness / (1 - nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    
    a1 = x2 * y3 - x3 * y2;
    a2 = x3 * y1 - x1 * y3;
    a3 = x1 * y2 - x2 * y1;
    
    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;
    
    c1 = x3 - x2;
    c2 = x1 - x3;
    c3 = x2 - x1;
    
    B = 1 / (2 * area) * [b1, 0, b2, 0, b3, 0; 0, c1, 0, c2, 0, c3; c1, b1, c2, b2, c3, b3];
    
    elementK = area * (B' * D * B);
end

function [globalK] = assemble_element_stiffness(globalK, elementK, elementNodes)
    % 组装单元刚度矩阵到全局刚度矩阵
    for i = 1:3
        for j = 1:3
            row = 2 * elementNodes(i) - 1;
            col = 2 * elementNodes(j) - 1;
            globalK(row:row+1, col:col+1) = globalK(row:row+1, col:col+1) + elementK(2*i-1:2*i, 2*j-1:2*j);
        end
    end
end