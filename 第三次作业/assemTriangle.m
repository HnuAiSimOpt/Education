%% S230200171郑小露 定义刚度矩阵组装
function k_t = assemTriangle(k_t, k_ele, node1, node2, node3)
    % assemTriangle 函数将平面三角形元素的单元刚度矩阵 k_ele 组装到全局刚度矩阵 k_t 中
    % 元素有节点 node1、node2 和 node3
    % 函数在组装后返回全局刚度矩阵 k_t

    % 计算节点自由度
    d(1:2) = 2 * node1 - 1:2 * node1;
    d(3:4) = 2 * node2 - 1:2 * node2;
    d(5:6) = 2 * node3 - 1:2 * node3;

    % 遍历局部刚度矩阵的所有元素，将其组装到全局刚度矩阵中
    for ii = 1:6
        for jj = 1:6
            k_t(d(ii), d(jj)) = k_t(d(ii), d(jj)) + k_ele(ii, jj);
        end
    end
end