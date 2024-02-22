

function [NDxyz,JacobiDET] = ShapeFunction(ElementNodeCoordinate)
%�����κ������κ����Ծֲ�����ksi eta zeta�ĵ���
NDL = [-1 1 0 0;-1 0 1 0;-1 0 0 1];%3*4  [N1Dksi N2Dksi N3Dksi N4Dksi;N1Deta N2Deta N3Deta N4Deta����]
Jacobi = NDL*ElementNodeCoordinate;%�����ſɱȾ���3*4 4*3
JacobiDET = det(Jacobi);%�����ſɱ�����ʽ3*3 [DxDksi DyDksi DzDksi;DxDeta����
JacobiINV=inv(Jacobi);%���ſɱ�����ʽ����3*3
NDxyz=JacobiINV*NDL;%�����ſɱ�����ʽ��������κ����Խṹ����ĵ���[DN1Dx DN2Dx DN3Dx;DN1Dy DN2Dy DN3Dy;����]
end