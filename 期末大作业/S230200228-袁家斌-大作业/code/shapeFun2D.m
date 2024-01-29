function [shapeMatrix, naturalDerivatives] = shapeFun2D(xi, eta, elem_type)

%% 该函数用于生成形函数矩阵和自然坐标系微分矩阵
% shapeMatrix 形函数矩阵 4 * 1(Q4) ; 3 * 1(T3) 
% naturalDerivatives 自然坐标导数矩阵 4 * 2(Q4) ; 3 * 2(T3) 
% xi, eta 自然坐标系下的坐标
% elem_type 单元类型

%% 
switch elem_type
    % 四节点四边形单元
    case 'Q4'
       % 形函数矩阵
        shapeMatrix = 0.25 * [(1-xi)*(1-eta); (1+xi)*(1-eta); ...
        (1+xi)*(1+eta); (1-xi)*(1+eta)];
       % 自然坐标导数矩阵
        naturalDerivatives = 0.25 * [-(1-eta), -(1-xi)
                            1-eta, -(1+xi)
                            1+eta, 1+xi
                            -(1+eta), 1-xi];
    % 三节点三角形单元
    case 'T3'
        shapeMatrix = [1-xi-eta; xi; eta];
        naturalDerivatives = [-1 -1; 1 0; 0 1];
end                   
end