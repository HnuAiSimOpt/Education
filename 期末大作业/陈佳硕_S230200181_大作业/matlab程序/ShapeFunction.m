%%%�����κ������κ��������ȫ������ĵ������ſɱ�����ʽ
%  N�ֲ������µ��κ�������   
%  NDerivative�κ��������ȫ������ĵ���  
%  JacobiDET�ſɱ�����ʽ
%  GaussPoint��˹������  
%  ElementNodeCoordinate��Ԫ�ڵ����꣨8*3��ÿһ�д���һ���ڵ�����꣩  
function [N,NDerivative,JacobiDET] = ShapeFunction(GaussPoint,ElementNodeCoordinate)
 %�Ȳ�Ԫ����  ÿһ�д���һ���������
ParentNodes=[-1  1  1 -1 -1  1  1 -1;
    -1 -1  1  1 -1 -1  1  1;
    -1 -1 -1 -1  1  1  1  1];
N=zeros(8,1); %��ʼ���κ�������8*1
ParentNDerivative=zeros(3,8);%��ʼ���κ����Ծֲ�����ĵ�������3*8
%�����κ������κ����Ծֲ����굼��
for I=1:8
    XPoint = ParentNodes(1,I);
    YPoint = ParentNodes(2,I);
    ZPoint = ParentNodes(3,I);
    ShapePart = [1+GaussPoint(1)*XPoint 1+GaussPoint(2)*YPoint 1+GaussPoint(3)*ZPoint];  %�����м����
    N(I) = 0.125*ShapePart(1)*ShapePart(2)*ShapePart(3);
    ParentNDerivative(1,I) = 0.125*XPoint*ShapePart(2)*ShapePart(3);
    ParentNDerivative(2,I) = 0.125*YPoint*ShapePart(1)*ShapePart(3);
    ParentNDerivative(3,I) = 0.125*ZPoint*ShapePart(1)*ShapePart(2);
end
Jacobi = ParentNDerivative*ElementNodeCoordinate;%�����ſɱȾ���
JacobiDET = det(Jacobi);%�����ſɱ�����ʽ
JacobiINV=inv(Jacobi);%���ſɱ�����ʽ����
NDerivative=JacobiINV*ParentNDerivative;%�����ſɱ�����ʽ��������κ����Խṹ����ĵ���
end