function [elemStiffness] = ...
    elemStiffness2D(elemcoords, thickness, constMatrix)

%% 此函数构成单元的刚度矩阵
% elemStiffness 单元刚度矩阵 24 * 24
% elemcoords 坐标矩阵 4 * 2 
% thickness 单元厚度值
% constMatrix 本构矩阵 3 * 3

%% 初始化
dof_m = [1 2 7 8 13 14 19 20];      % 膜效应对应自由度
dof_b = [4 5 10 11 16 17 22 23];    % 弯曲对应自由度
dof_s = [3 4 5 9 10 11 15 16 17 21 22 23];  % 横向剪切对应自由度
elemStiffness = zeros(24);

%% 计算并整合膜效应和平面弯曲刚度矩阵
[weights, points] = gaussIntegration(2);
for g = 1:length(weights)
    xi = points(g, 1);
    eta = points(g, 2);
    [~, dN] = shapeFun2D(xi, eta, 'Q4');
    [J, dXY] = jacobian2D(elemcoords, dN);
    [Bm, Bb] = strainMatrix2D(dXY);
    J_delta = det(J);
    elemStiffness(dof_m, dof_m) = elemStiffness(dof_m, dof_m)...
    + thickness * Bm' * constMatrix.in * Bm * J_delta * weights(g);

     elemStiffness(dof_b, dof_b) = elemStiffness(dof_b, dof_b)...
    + thickness^3 * Bb' * constMatrix.in * Bb * J_delta * weights(g) / 12;
end

%% 计算并整合横向剪切刚度矩阵
[weights, points] = gaussIntegration(1);
for g = 1:length(weights)
    xi = points(g, 1);
    eta = points(g, 2);
    [N, dN] = shapeFun2D(xi, eta, 'Q4');
    [J, dXY] = jacobian2D(elemcoords, dN);
    [~, ~, Bs] = strainMatrix2D(dXY, N);
    
    elemStiffness(dof_s, dof_s) = elemStiffness(dof_s, dof_s) + ...
    constMatrix.alpha * thickness * Bs' * constMatrix.out * Bs * det(J) * weights(g);
end

end