function [JacobianMatrix, PysicalDerivatives] = ...
    jacobian2D(elemCoordinates, naturalDerivatives)

%% �˺������������ſɱȾ�����κ�������������ĵ�������
% JacobianMatrix �ſɱȾ��� 2 * 2
% PysicalDerivatives�κ�������������ĵ������� 2 * 4(Q4); 2 * 3(T3)
% elemCoordinates ��Ԫ�������4 * 2(Q4); 3 * 2(T3)
% naturalDerivatives ��Ȼ���굼������ ��Ϊ 4 * 2(Q4); 3 * 2(T3)

JacobianMatrix = naturalDerivatives' * elemCoordinates;
PysicalDerivatives = JacobianMatrix \ naturalDerivatives';
end

