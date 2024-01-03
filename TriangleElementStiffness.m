%% S230200171郑小露 定义单元刚度矩阵
function k_ele=TriangleElementStiffness(E,v,t,node_ele)
    % S230200171 郑小露
    % 定义元素刚度矩阵。
    % 元素刚度矩阵的大小是6 x 6。

    %-----------% 节点坐标 %--------%
    x1 = node_ele(1,1);                
    y1 = node_ele(1,2);
    x2 = node_ele(2,1);                
    y2 = node_ele(2,2);
    x3 = node_ele(3,1);                
    y3 = node_ele(3,2);
    %-------------------------------%

    A = 1/2*[1 x1 y1;1 x2 y2;1 x3 y3]
    A = det(A);   % 单元面积
    a1 = [x2 y2;x3 y3];a1 = det(a1);
    a2 = -[x1 y1;x3 y3];a2 = det(a2);
    a3 = [x1 y1;x2 y2];a3 = det(a3);
    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;
    c1 = x3 - x2;
    c2 = x1 - x3;
    c3 = x2 - x1;
    B = 1/2/A * [b1 0 b2 0 b3 0;
             0 c1 0 c2 0 c3;
             c1 b1 c2 b2 c3 b3];
    D = E/(1-v^2) * [1 v 0;
                        v 1 0;
                        0 0 (1-v)/2];     % 平面应力矩阵

    k_ele=B'*D*B*t*A;    % 单元刚度矩阵
end