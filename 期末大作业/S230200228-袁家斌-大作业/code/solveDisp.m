function [displacement] = solveDisp(stiffness, disp, load)

%% �˺�����֪�غ���λ��
% displacement �ڵ��������ɶ�λ��
% stiffness ȫ�ָնȾ���
% disp Լ�������ɶ�
% load �ڵ��غ�

%% 
stiffness(disp, :) = 0;
stiffness(:, disp) = 0;
stiffness(disp, disp) = eye(length(disp));
load(disp) = 0;
displacement = stiffness \ load;
end