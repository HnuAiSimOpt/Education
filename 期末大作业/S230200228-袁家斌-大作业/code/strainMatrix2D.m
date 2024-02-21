function [Bm, Bb, Bs] = strainMatrix2D(PysicalDerivatives, shapeFun)

%% �˺�����������Ӧ�����
% Bm ����ĤӦ����� 3 * 8
% Bb ��������Ӧ����� 3 * 8
% Bs �������Ӧ����� 3 * 12
% PysicalDerivatives �κ������������굼������ 2 * 4
% shapeFun �κ�������

if nargin == 1
    %% ����ĤӦ�����
    Bm = zeros(3, 8);
    Bm(1, 1:2:end) = PysicalDerivatives(1,:);
    Bm(2, 2:2:end) = PysicalDerivatives(2,:);
    Bm(3, 1:2:end) = PysicalDerivatives(2,:);
    Bm(3, 2:2:end) = PysicalDerivatives(1,:);

    %% ��������Ӧ�����
    Bb = zeros(3, 8);
    Bb(1, 2:2:end) = PysicalDerivatives(1,:);
    Bb(2, 1:2:end) = -PysicalDerivatives(2,:);
    Bb(3, 1:2:end) = -PysicalDerivatives(1,:);
    Bb(3, 2:2:end) = PysicalDerivatives(2,:);
    
    Bs = 0;
else
    %% �������Ӧ�����
    Bs = zeros(2, 12);
    Bs(1, 1:3:end) = PysicalDerivatives(1,:);
    Bs(1, 3:3:end) = shapeFun;
    Bs(2, 1:3:end) = PysicalDerivatives(2,:);
    Bs(2, 2:3:end) = -shapeFun;
    
    Bm = 0; Bb = 0;
end

end