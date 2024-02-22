% 计算应变场
epsilon_x = zeros(size(elements, 1), 1);
epsilon_y = zeros(size(elements, 1), 1);

for e = 1:size(elements, 1)
    nodes_e = elements(e, :);
    x_e = nodes(nodes_e, 1);
    y_e = nodes(nodes_e, 2);

    % 计算单元面积
    A_e = 0.5 * abs(det([1, 1, 1; x_e(:)'; y_e(:)']));

    % 计算 b_i 和 c_i
    b = [y_e(2) - y_e(3), y_e(3) - y_e(1), y_e(1) - y_e(2)];
    c = [x_e(3) - x_e(2), x_e(1) - x_e(3), x_e(2) - x_e(1)];

    % 计算 B 矩阵
    B = 1 / (2 * A_e) * [b(1), 0, b(2), 0, b(3), 0;
                        0, c(1), 0, c(2), 0, c(3);
                        c(1), b(1), c(2), b(2), c(3), b(3)];

    % 计算单元位移
    u_e = displacements(2*nodes_e - 1);
    v_e = displacements(2*nodes_e);

    % 计算应变场
    epsilon_x(e) = B(1, :) * [u_e; v_e];
    epsilon_y(e) = B(2, :) * [u_e; v_e];
end

% 计算应力场
sigma_x = E / (1 - nu^2) * (epsilon_x + nu * epsilon_y);
sigma_y = E / (1 - nu^2) * (nu * epsilon_x + epsilon_y);

% 绘制应力场颜色图
figure;
pdeplot(mesh, 'XYData', sigma_x, 'ColorMap', 'jet', 'XYStyle', 'interp');
title('X方向应力场');
xlabel('X轴位置');
ylabel('Y轴位置');
colorbar;

figure;
pdeplot(mesh, 'XYData', sigma_y, 'ColorMap', 'jet', 'XYStyle', 'interp');
title('Y方向应力场');
xlabel('X轴位置');
ylabel('Y轴位置');
colorbar;
