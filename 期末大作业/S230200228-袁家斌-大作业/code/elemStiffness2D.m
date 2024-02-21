function [elemStiffness] = ...
    elemStiffness2D(elemcoords, thickness, constMatrix)

%% �˺������ɵ�Ԫ�ĸնȾ���
% elemStiffness ��Ԫ�նȾ��� 24 * 24
% elemcoords ������� 4 * 2 
% thickness ��Ԫ���ֵ
% constMatrix �������� 3 * 3

%% ��ʼ��
dof_m = [1 2 7 8 13 14 19 20];      % ĤЧӦ��Ӧ���ɶ�
dof_b = [4 5 10 11 16 17 22 23];    % ������Ӧ���ɶ�
dof_s = [3 4 5 9 10 11 15 16 17 21 22 23];  % ������ж�Ӧ���ɶ�
elemStiffness = zeros(24);

%% ���㲢����ĤЧӦ��ƽ�������նȾ���
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

%% ���㲢���Ϻ�����иնȾ���
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