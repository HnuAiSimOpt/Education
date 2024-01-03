function k_t = assemTriangle(k_t, k_ele, node1 ,node2, node3)
    % assemTriangle 这个函数将平面三角形元素的元素刚度矩阵k组装到全局刚度矩阵K中。
    % 元素有节点i和j。
    % 函数在元素刚度矩阵k被组装后返回全局刚度矩阵K。

    d(1:2) = 2*node1 - 1:2*node1;
    d(3:4) = 2*node2 - 1:2*node2;
    d(5:6) = 2*node3 - 1:2*node3;
    for ii = 1 : 6
        for jj = 1 : 6
            k_t(d(ii), d(jj)) = k_t(d(ii), d(jj)) + k_ele(ii, jj);
        end
    end
end