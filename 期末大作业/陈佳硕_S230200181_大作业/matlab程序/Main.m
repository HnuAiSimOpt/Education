%%%������
%%%������Ԫ��������T������������ѧ������һ�������嵥Ԫ��
%%%�������¼�˶
%%%ѧ�ţ�S230200181
clc;clear;
%��ȡ�ļ���ýڵ�������ϢNodes����Ԫ��ϢElements
[Nodes, Elements] = Readmessage( 'T.inp' );
Nodes=Nodes*0.001;%��λ���㣨mm--m��
Mg=5000;%����T���������غɣ�N��
% �������� Forces=[�����ڵ�  ��������(1,2,3�ֱ����x,y,z)  ������С]
row_number=6;%��������Ľڵ�����
arrange_number=12;%��������Ľڵ�����
gap=77;%ÿ�м������
Forces=zeros(row_number*arrange_number,3);
Forces(:,2)=3;%��������Ϊz�᷽��
Forces(:,3)=-(Mg/(row_number*arrange_number));%ÿ�������ڵ���غɣ����ű�ʾz�Ḻ����
for i=1:row_number
    Forces(1+(i-1)*arrange_number,1)=1+(i-1)*gap;
    Forces(2+(i-1)*arrange_number,1)=14+(i-1)*gap;
    Forces(3+(i-1)*arrange_number:3+(i-1)*arrange_number+9,1)=41+(i-1)*gap:41+(i-1)*gap+9;
end
%Լ���ڵ�ı��
ConNumber=3851:1:3927;
%Լ������ Constraints=[ǿ��λ�ƽڵ�  ǿ��λ�Ʒ���(1,2,3�ֱ����x,y,z)  ǿ��λ�ƴ�С]  
Constraints=zeros(size(ConNumber,2)*3,3);
for i=1:size(ConNumber,2)
Constraints(3*i-2:3*i,:)=[ConNumber(i) 1 0;ConNumber(i) 2 0;ConNumber(i) 3 0;];
end
E=2.1e11; %����ģ��
u=0.3;    %���ɱ�
%���ú������Ӧ��Ӧ������D
D=LinearIsotropicD(E,u);
%���λ��U
Dof=3;
NodeCount = size(Nodes,1);  % �ڵ����
ElementCount= size(Elements,1);  %��Ԫ����
Dofs = Dof*NodeCount; %�����ɶ���
U=sparse(Dofs,1);	% ��ʼ���ṹλ��
K = sparse(Dofs,Dofs);  %��ʼ������ն���
Force = sparse(Dofs,1);  %��ʼ����������
%����Ӧ��-Ӧ�����D
D=LinearIsotropicD(E,u);
for I=1:ElementCount
    % ��Ԫ�ڵ�����
    ElementNodeCoordinate=Nodes(Elements(I,:),:);
    % ���㵥��
    ElementStiffnessMatrix=Ke(D,ElementNodeCoordinate);
    % ���㵥Ԫ�ڵ����ɶȱ��
    ElementNodeDOF=zeros(1,24);
    for J=1:8
        II=(J-1)*Dof+1;
        ElementNodeDOF(II:II+2)=(Elements(I,J)-1)*Dof+1:(Elements(I,J)-1)*Dof+3;
    end
    K(ElementNodeDOF,ElementNodeDOF)=K(ElementNodeDOF,ElementNodeDOF)+ElementStiffnessMatrix;%���ݵ�Ԫ�ڵ����ɶ�  ��װ�ܸ�
end
% ʩ������
if size(Forces,1)>0
    ForceDOF = Dof*(Forces(:,1)-1)+Forces(:,2);   %�����������ɶȱ��
    Force(ForceDOF) = Force(ForceDOF) + Forces(:,3);
end
% �˴�����ʩ��λ��Լ��
BigNumber=1e8;
ConstraintsNumber=size(Constraints,1);
if ConstraintsNumber~=0
    FixedDof=Dof*(Constraints(:,1)-1)+Constraints(:,2);  %��Լ�������ɶȱ��(������)
    for i=1:ConstraintsNumber
        K(FixedDof(i),FixedDof(i))=K(FixedDof(i),FixedDof(i))*BigNumber;
        Force(FixedDof(i))=Constraints(i,3)*K(FixedDof(i),FixedDof(i));
    end
end
%����λ��
U = K\Force;
% ���Ӧ����Ӧ�䲢�����Ӧ��ͼ
OutputResults(Nodes,Elements,D,U)%�����������ļ�