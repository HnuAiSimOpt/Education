%% S230200171郑小露 定义单元应力矩阵
function str = TriangleElementStress(E, v, node_ele, u1, p)
    % TriangleElementStress 函数返回具有弹性模量 E、泊松比 v、
    % 元素的节点坐标 node_ele 的三角形元素的元素应力矩阵，平面应力或平面应变选项。

    %-----------% 节点坐标 %--------%
    x1 = node_ele(1, 1);
    y1 = node_ele(1, 2);
    x2 = node_ele(2, 1);
    y2 = node_ele(2, 2);
    x3 = node_ele(3, 1);
    y3 = node_ele(3, 2);
    %-------------------------------%

    % 计算三角形面积
    A = 1/2 * [1 x1 y1; 1 x2 y2; 1 x3 y3];
    A = det(A);

    % 计算向量 a1, a2, a3
    a1 = det([x2 y2; x3 y3]);
    a2 = -det([x1 y1; x3 y3]);
    a3 = det([x1 y1; x2 y2]);

    % 计算向量 b1, b2, b3
    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;

    % 计算向量 c1, c2, c3
    c1 = x3 - x2;
    c2 = x1 - x3;
    c3 = x2 - x1;

    % 计算形函数矩阵 B
    B = 1/2/A * [b1 0 b2 0 b3 0;
                 0 c1 0 c2 0 c3;
                 c1 b1 c2 b2 c3 b3];

    % 根据平面应力或平面应变选项选择弹性矩阵 D
    if p == 1
        D = E / (1 - v^2) * [1 v 0;
                             v 1 0;
                             0 0 (1 - v)/2]; % 平面应力矩阵
    elseif p == 2
        D = E / (1 + v) / (1 - 2 * v) * [1 - v v 0;
                                          v 1 - v 0;
                                          0 0 (1 - 2 * v)/2]; % 平面应变矩阵
    end
    size_D = size(D);
    size_B = size(B);
    size_u1 = size(u1);
    disp(['Size of D: ', num2str(size_D)]);
    disp(['Size of B: ', num2str(size_B)]);
    disp(['Size of u1: ', num2str(size_u1)]);
    % 计算单元应力
    str = D * B * u1;
end
