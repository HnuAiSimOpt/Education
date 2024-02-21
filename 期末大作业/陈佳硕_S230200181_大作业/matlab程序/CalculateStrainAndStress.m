%%%����Ӧ��Ӧ��(�ڵ㼰��˹��)
%  NodeStrain�ڵ�Ӧ��
%  NodeStress�ڵ�Ӧ��  
%  GaussStrain��˹��Ӧ��   
%  GaussStress��˹��Ӧ��
%  Uλ�ƾ���
%  D�ߵ��Բ���Ӧ��-Ӧ�����
%  Nodes�ڵ�������Ϣ
%  Elements��Ԫ��Ϣ
function [NodeStrain,NodeStress,GaussStrain,GaussStress]=CalculateStrainAndStress(U,D,Nodes,Elements)  
ElementCount= size(Elements,1);  %��Ԫ����
GaussCoordinate=[-0.57735026918963D0, 0.57735026918963D0];  %��˹���ֵ�����
GaussWeight=[1.00000000000000D0, 1.00000000000000D0];  %��˹���ֵ�Ȩ��
GaussPointNumber=0;%��˹���ֵ���
INODE=0;% �ڵ���
Dof=3; 
%�������ά�� �ӿ������ٶ�
GaussStrain=zeros(6,ElementCount*8);
GaussStress=zeros(6,ElementCount*8);
NodeStrain=zeros(6,ElementCount*8);
NodeStress=zeros(6,ElementCount*8);
%ѭ����װ�ܸ�
for I=1:ElementCount
    % ��Ԫ�ڵ�����
    ElementNodeCoordinate=Nodes(Elements(I,:),:); 
    % ���㵥Ԫ�ڵ����ɶȱ��
    ElementNodeDOF=zeros(1,24);
    for J=1:8
        II=(J-1)*Dof+1;
        ElementNodeDOF(II:II+2)=(Elements(I,J)-1)*Dof+1:(Elements(I,J)-1)*Dof+3;
    end
    K=1;
    InterpolationMatrix=zeros(8,8);%���ڵ�Ӧ��Ӧ��Ĳ�ֵ����
    %ѭ����˹��
    for X=1:2
        for Y=1:2
            for Z=1:2  
                E1=GaussCoordinate(X); E2=GaussCoordinate(Y); E3=GaussCoordinate(Z);
                GaussPointNumber = GaussPointNumber + 1;
                % ����ֲ������µ��κ������κ�������
                [N,NDerivative, ~] = ShapeFunction([E1 E2 E3], ElementNodeCoordinate);
                ElementNodeDisplacement=U(ElementNodeDOF);
                ElementNodeDisplacement=reshape(ElementNodeDisplacement,Dof,8);
                % �����˹��Ӧ�� GausspointStrain3_3 3*3��Ӧ�����
                GausspointStrain3_3=ElementNodeDisplacement*NDerivative';
                %�Ѹ�˹��Ӧ������д��6*1
                GausspointStrain=[GausspointStrain3_3(1,1) GausspointStrain3_3(2,2) GausspointStrain3_3(3,3) ...
                GausspointStrain3_3(1,2)+GausspointStrain3_3(2,1)....
                GausspointStrain3_3(2,3)+GausspointStrain3_3(3,2) GausspointStrain3_3(1,3)+GausspointStrain3_3(3,1)]';
                % �����˹��Ӧ��
                GausspointStress = D*GausspointStrain;
                GaussStrain(:,GaussPointNumber)=GausspointStrain;
                GaussStress(:,GaussPointNumber)=GausspointStress;
                InterpolationMatrix(K,:)=N;
                K=K+1;
            end
        end
    end
    %��ýڵ�Ӧ��Ӧ��
    Temp1=InterpolationMatrix\(GaussStrain(1:6,GaussPointNumber-7:GaussPointNumber)');
    NodeStrain(1:6,INODE+1:INODE+8)=Temp1';
    Temp2=InterpolationMatrix\(GaussStress(1:6,GaussPointNumber-7:GaussPointNumber)');
    NodeStress(1:6,INODE+1:INODE+8)=Temp2';
    INODE=INODE+8;
end
end