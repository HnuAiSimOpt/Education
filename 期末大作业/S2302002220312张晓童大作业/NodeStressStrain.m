
function [NodeStrain,NodeStress]=NodeStressStrain(U,D,Nodes,Elements)  
ElementNum= size(Elements,1);  %��Ԫ����
Dof=3; 
NodeStrain=zeros(6,ElementNum*4);%11 22 33 12 13 23
NodeStress=zeros(6,ElementNum*4);
Strain = zeros(6,1);
Stress = zeros(6,1);
%ѭ����װ�ܸ�
for i=1:ElementNum
    % ��Ԫ�ڵ�����
    ElementNodeCoordinate=Nodes(Elements(i,:),:); %4*3
    % ���㵥Ԫ�ڵ����ɶȱ��,������λ��
    ElementNodeDOF=zeros(1,12);%3*4
    for j=1:4
        ii=(j-1)*Dof+1;
        %��Ԫ���ڵ����ܸնȾ����е�λ������
        ElementNodeDOF(ii:ii+2)=(Elements(i,j)-1)*Dof+1:(Elements(i,j)-1)*Dof+3;
    end
    % �����κ�������
    [NDxyz, ~] = ShapeFunction(ElementNodeCoordinate);%[DN1Dx DN2Dx DN3Dx;DN1Dy DN2Dy DN3Dy;����]
    ElementNodeDisplacement=U(ElementNodeDOF);%12*1 �ڵ�λ������
    ElementNodeDisplacement=reshape(ElementNodeDisplacement,Dof,4);%3*4
    % ���㵥ԪӦ�� Strain3_3 3*3��Ӧ�����
    Strain3_3=ElementNodeDisplacement*NDxyz';%��˹���ֵ㴦Ӧ�� 3*4  4*3
    %�ѵ�ԪӦ������д��6*1
    Strain=[Strain3_3(1,1) Strain3_3(2,2) Strain3_3(3,3) ...
    Strain3_3(1,2)+Strain3_3(2,1)....
    Strain3_3(2,3)+Strain3_3(3,2) Strain3_3(1,3)+Strain3_3(3,1)]';
    Stress(1:6,1) = D*Strain;%��˹���ֵ㴦Ӧ��
    %��ýڵ�Ӧ��Ӧ��
    for X = 1:4
        NodeStrain(1:6,((i-1)*4+X)) = Strain(1:6,1);%���ո���Ԫ4�ڵ���������1�ŵ�Ԫ4�ڵ㣻2�ŵ�Ԫ4�ڵ�
        NodeStress(1:6,((i-1)*4+X)) = Stress(1:6,1);
    end
end
end