function [JacobianMatrix, PysicalDerivatives] = ...
    jacobian2D(elemCoordinates, naturalDerivatives)

%% 此函数用于生成雅可比矩阵和形函数对物理坐标的导数矩阵
% JacobianMatrix 雅可比矩阵 2 * 2
% PysicalDerivatives形函数对物理坐标的导数矩阵 2 * 4(Q4); 2 * 3(T3)
% elemCoordinates 单元坐标矩阵，4 * 2(Q4); 3 * 2(T3)
% naturalDerivatives 自然坐标导数矩阵 ，为 4 * 2(Q4); 3 * 2(T3)

JacobianMatrix = naturalDerivatives' * elemCoordinates;
PysicalDerivatives = JacobianMatrix \ naturalDerivatives';
end

