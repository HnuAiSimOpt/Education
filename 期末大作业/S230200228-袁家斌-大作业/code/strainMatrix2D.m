function [Bm, Bb, Bs] = strainMatrix2D(PysicalDerivatives, shapeFun)

%% 此函数用于生产应变矩阵
% Bm 面内膜应变矩阵 3 * 8
% Bb 面外弯曲应变矩阵 3 * 8
% Bs 面外剪切应变矩阵 3 * 12
% PysicalDerivatives 形函数对物理坐标导数矩阵 2 * 4
% shapeFun 形函数矩阵

if nargin == 1
    %% 面内膜应变矩阵
    Bm = zeros(3, 8);
    Bm(1, 1:2:end) = PysicalDerivatives(1,:);
    Bm(2, 2:2:end) = PysicalDerivatives(2,:);
    Bm(3, 1:2:end) = PysicalDerivatives(2,:);
    Bm(3, 2:2:end) = PysicalDerivatives(1,:);

    %% 面外弯曲应变矩阵
    Bb = zeros(3, 8);
    Bb(1, 2:2:end) = PysicalDerivatives(1,:);
    Bb(2, 1:2:end) = -PysicalDerivatives(2,:);
    Bb(3, 1:2:end) = -PysicalDerivatives(1,:);
    Bb(3, 2:2:end) = PysicalDerivatives(2,:);
    
    Bs = 0;
else
    %% 面外剪切应变矩阵
    Bs = zeros(2, 12);
    Bs(1, 1:3:end) = PysicalDerivatives(1,:);
    Bs(1, 3:3:end) = shapeFun;
    Bs(2, 1:3:end) = PysicalDerivatives(2,:);
    Bs(2, 2:3:end) = -shapeFun;
    
    Bm = 0; Bb = 0;
end

end