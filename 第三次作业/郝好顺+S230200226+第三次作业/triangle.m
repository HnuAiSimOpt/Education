%郝好顺，S230200226
%采用两个三角形网格计算一矩形薄板的应力，边界条件左端固定，右端水平方向力为0

clear;
clc;

%参数
a = 4;
b = 3;
t = 10e-3;
E = 210e9;
mu = 0.1;
F = 8e6;
h = 0.2;

%划分网格，生成节点和单元信息
[Node_info, Ele_info] = Meshing(a, b, h);
%计算单元刚度矩阵，组装总刚度矩阵
[K, D, BB] = Assembly(Node_info, Ele_info, E, mu, t);
%生成节点等效载荷
R = Load(Node_info, F);
%引入边界条件
[KK, RR] = BC(Node_info, K, R);
u = lsqminnorm(KK, RR);
%计算每个单元的应力
Stress(BB, Ele_info, D, u);

function [Node_info, Ele_info] = Meshing(a, b, h)
    Num_nodes = 0;
    Node_info = [];

    for i = 0:h:a
        for j = 0:h:b
            Num_nodes = Num_nodes + 1;
            Node_info = [Node_info; Num_nodes, i, j];
        end
    end

    Num_eles = 0;
    Ele_info = [];

    for i = h:h:a
        for j = h:h:b
            Num_eles = Num_eles + 1;
            node_local = [i-h, j-h; i, j-h; i-h, j];
            node_list = [];

            for k = 1:3
                row_x = find(abs(Node_info(:,2) - node_local(k, 1)) < 1e-6);
                row_y = find(abs(Node_info(:,3) - node_local(k, 2)) < 1e-6);
                num_node = intersect(row_x, row_y);
                node_list = [node_list; num_node];
            end

            Ele_info = [Ele_info; Num_eles, node_list'];
            Num_eles = Num_eles + 1;
            node_local = [i, j; i-h, j; i, j-h];
            node_list = [];

            for k = 1:3
                row_x = find(abs(Node_info(:,2) - node_local(k, 1)) < 1e-6);
                row_y = find(abs(Node_info(:,3) - node_local(k, 2)) < 1e-6);
                num_node = intersect(row_x, row_y);
                node_list = [node_list; num_node];
            end

            Ele_info = [Ele_info; Num_eles, node_list'];
        end
    end
end

function [K, D, BB] = Assembly(Node_info, Ele_info, E, mu, t)
    Num_nodes = size(Node_info, 1);
    Num_eles = size(Ele_info, 1);
    K = zeros(2 * Num_nodes);
    BB = [];

    for i = 1:Num_eles
        node_info_local = [Ele_info(i, 2), Node_info(Ele_info(i, 2), 2:3); ...
                           Ele_info(i, 3), Node_info(Ele_info(i, 3), 2:3); ...
                           Ele_info(i, 4), Node_info(Ele_info(i, 4), 2:3)];
        [ke, D, B] = Ke(node_info_local, E, mu, t);
        BB = [BB; B]; 

        j = node_info_local(1, 1);
        k = node_info_local(2, 1);
        m = node_info_local(3, 1);
        num = [2 * j-1, 2 * j, 2 * k-1, 2 * k, 2 * m-1, 2 * m];

        for n1 = 1:6
            for n2 = 1:6
                K(num(n1), num(n2)) = K(num(n1), num(n2)) + ke(n1, n2);
            end
        end
    end
end

function [ke, D, B] = Ke(node_info, E, mu, t)
    C = [1, node_info(1, 2), node_info(1, 3); ...
         1, node_info(2, 2), node_info(2, 3); ...
         1, node_info(3, 2), node_info(3, 3)];
    A = 0.5 * det(C);
    B = 0.5 / A * [node_info(2, 3) - node_info(3, 3), 0, node_info(3, 3) - node_info(1, 3), 0, node_info(1, 3) - node_info(2, 3), 0; ...
                   0, node_info(3, 2) - node_info(2, 2), 0, node_info(1, 2) - node_info(3, 2), 0, node_info(2, 2) - node_info(1, 2); ...
                   node_info(3, 2) - node_info(2, 2), node_info(2, 3) - node_info(3, 3), node_info(1, 2) - node_info(3, 2), ...
                   node_info(3, 3) - node_info(1, 3), node_info(2, 2) - node_info(1, 2), node_info(1, 3) - node_info(2, 3)];
    D = E / (1 - mu^2) * [1, mu, 0; mu, 1, 0; 0, 0, (1 - mu) / 2];
    ke = B' * D * B * A * t;
end

function R = Load(Node_info, F)
    R = zeros(size(Node_info, 1) * 2, 1);
    row1_x = find(Node_info(:, 2) == 2);
    row1_y = find(Node_info(:, 3) == 0);
    num1 = intersect(row1_x, row1_y);
    row2_x = find(Node_info(:, 2) == 2);
    row2_y = find(Node_info(:, 3) == 1);
    num2 = intersect(row2_x, row2_y);
    R(2 * num1) = -F / 2;
    R(2 * num2) = -F / 2;
end

function [KK, RR] = BC(Node_info, K, R)
    num = find(Node_info(:, 2) < 1e-9);
    KK = K;
    RR = R;

    for i = 1:size(num, 1)
        r = num(i);
        KK(2 * r - 1, :) = 0;
        KK(:, 2 * r - 1) = 0;
        KK(2 * r - 1, 2 * r - 1) = 1;
        KK(2 * r, :) = 0;
        KK(:, 2 * r) = 0;
        KK(2 * r, 2 * r) = 1;
        RR(2 * r - 1) = 0;
        RR(2 * r) = 0;
    end
end

function Stress(BB, Ele_info, D, u)
    Num_eles = size(Ele_info, 1);
    sigma_x = zeros(2, 1);
    sigma_y = zeros(2, 1);
    sigma_xy = zeros(2, 1);

    % 最右上角节点
    node_local = Ele_info(end-1, 2:4);
    u_local = [u(2 * node_local(1) - 1); u(2 * node_local(1)); ...
               u(2 * node_local(2) - 1); u(2 * node_local(2)); ...
               u(2 * node_local(3) - 1); u(2 * node_local(3))];
    stress_element = D * BB(3 * (Num_eles-1) - 2:3 * (Num_eles-1), :) * u_local;
    sigma_x(1) = stress_element(1) / 1e9;
    sigma_y(1) = stress_element(2) / 1e9;
    sigma_xy(1) = stress_element(3) / 1e9;

    % 最右下角节点
    node_local = Ele_info(end, 2:4);
    u_local = [u(2 * node_local(1) - 1); u(2 * node_local(1)); ...
               u(2 * node_local(2) - 1); u(2 * node_local(2)); ...
               u(2 * node_local(3) - 1); u(2 * node_local(3))];
    stress_element = D * BB(3 * Num_eles - 2:3 * Num_eles, :) * u_local;
    sigma_x(2) = stress_element(1) / 1e9;
    sigma_y(2) = stress_element(2) / 1e9;
    sigma_xy(2) = stress_element(3) / 1e9;

    disp('Stress Components for the right upper and lower nodes:');
    disp('Node | Sigma_x (GPa) | Sigma_y (GPa) | Sigma_xy (GPa)');
    disp([Ele_info(end-1:end, 1), sigma_x, sigma_y, sigma_xy]);
end

