function [displacement] = solveDisp(stiffness, disp, load)

%% 此函数已知载荷求位移
% displacement 节点所有自由度位移
% stiffness 全局刚度矩阵
% disp 约束的自由度
% load 节点载荷

%% 
stiffness(disp, :) = 0;
stiffness(:, disp) = 0;
stiffness(disp, disp) = eye(length(disp));
load(disp) = 0;
displacement = stiffness \ load;
end