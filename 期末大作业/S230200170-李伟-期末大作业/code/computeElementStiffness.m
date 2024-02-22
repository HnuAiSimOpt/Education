function Ke = computeElementStiffness(E, nu, t, x, y)
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
end
