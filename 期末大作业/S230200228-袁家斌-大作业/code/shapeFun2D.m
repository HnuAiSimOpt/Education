function [shapeMatrix, naturalDerivatives] = shapeFun2D(xi, eta, elem_type)

%% �ú������������κ����������Ȼ����ϵ΢�־���
% shapeMatrix �κ������� 4 * 1(Q4) ; 3 * 1(T3) 
% naturalDerivatives ��Ȼ���굼������ 4 * 2(Q4) ; 3 * 2(T3) 
% xi, eta ��Ȼ����ϵ�µ�����
% elem_type ��Ԫ����

%% 
switch elem_type
    % �Ľڵ��ı��ε�Ԫ
    case 'Q4'
       % �κ�������
        shapeMatrix = 0.25 * [(1-xi)*(1-eta); (1+xi)*(1-eta); ...
        (1+xi)*(1+eta); (1-xi)*(1+eta)];
       % ��Ȼ���굼������
        naturalDerivatives = 0.25 * [-(1-eta), -(1-xi)
                            1-eta, -(1+xi)
                            1+eta, 1+xi
                            -(1+eta), 1-xi];
    % ���ڵ������ε�Ԫ
    case 'T3'
        shapeMatrix = [1-xi-eta; xi; eta];
        naturalDerivatives = [-1 -1; 1 0; 0 1];
end                   
end