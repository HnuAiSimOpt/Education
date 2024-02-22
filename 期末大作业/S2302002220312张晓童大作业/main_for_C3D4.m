
clear all;
clc;
close all;
E=210000; %����ģ��
u=0.3;    %���ɱ�
D=E/((1+u)*(1-2*u))*[1-u u u 0 0 0;......
    u 1-u u 0 0 0;......
    u u 1-u 0 0 0;......
    0 0 0 (1-2*u)/2 0 0;......
    0 0 0 0 (1-2*u)/2 0;......
    0 0 0 0 0 (1-2*u)/2];%elastic modules matrix
load nodes_and_elements.mat
Forces=[7 2 -100;10 2 -100;12 2 -100;13 2 -100;14 2 -100;]; %[nodeID direction value]
ConID=[2,   3,   5,   8,  36,  37,  38,  61,  62,  63,  64,  87,  88,  89,  90,  91......
  92,  93, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260];
Constraints=zeros(size(ConID,2)*3,3);
for i=1:size(ConID,2)
Constraints(3*i-2:3*i,:)=[ConID(i) 1 0;ConID(i) 2 0;ConID(i) 3 0;];%[nodeID direction value];
end
Dof=3;                          % ÿ���ڵ��������ɶ�
NodeNum = size(Nodes,1);        % �ڵ����
ElementNum= size(Elements,1);   % ��Ԫ����
Dofs = Dof*NodeNum;             % �����ɶ���
U=sparse(Dofs,1);	            % ��ʼ���ṹλ�ƣ��վ��󣩣�����ϡ�������߼����ٶ�
K = sparse(Dofs,Dofs);          % ��ʼ������ն��󣨿վ���
F = sparse(Dofs,1);             % ��ʼ�������������վ���

for I=1:ElementNum
    % ��Ԫ�ڵ�����
    ElementNodeCoordinate=Nodes(Elements(I,:),:);
    % ���㵥Ԫ�նȾ���
    ElementStiffnessMatrix=Ke_C3D4(D,ElementNodeCoordinate);
    % ���㵥Ԫ�ڵ����ɶȱ��
    ElementNodeDOF=zeros(1,12);
    for J=1:4
        II=(J-1)*Dof+1;
        ElementNodeDOF(II:II+2)=(Elements(I,J)-1)*Dof+1:(Elements(I,J)-1)*Dof+3;
    end
    %����նȾ�����װ�����ݵ�Ԫ�ڵ����ɶ���װ�ܸ�
    K(ElementNodeDOF,ElementNodeDOF)=K(ElementNodeDOF,ElementNodeDOF)+ElementStiffnessMatrix;
end
% ʩ������
if size(Forces,1)>0
    %�����������ɶȱ��,��һ��Ϊ�ڵ�ţ��ڶ���Ϊ������
    ForceDOF = Dof*(Forces(:,1)-1)+Forces(:,2);
    F(ForceDOF) = F(ForceDOF) + Forces(:,3);
end
% �˴�����������������ʩ��λ��Լ��
BigNumber=1e8;
ConsNum=size(Constraints,1);
if ConsNum~=0
    FixedDof=Dof*(Constraints(:,1)-1)+Constraints(:,2);  %��Լ�������ɶȱ��(������)
    for i=1:ConsNum
        K(FixedDof(i),FixedDof(i))=K(FixedDof(i),FixedDof(i))*BigNumber;
        F(FixedDof(i))=Constraints(i,3)*K(FixedDof(i),FixedDof(i));
    end
end
% ����λ��
U = K\F;
% ������
% �����˹�㡢�ڵ�Ӧ��Ӧ��ֵ
[NodeStrain,NodeStress]=NodeStressStrain(U,D,Nodes,Elements);
% ���MISESӦ������
ElementNodeCount=4; % ÿ����Ԫ�ڵ���
MISES=zeros(1,ElementNum*ElementNodeCount);
for I=1:ElementNum*ElementNodeCount
    MISES(I)=sqrt(0.5)*sqrt((NodeStress(1,I)-NodeStress(2,I))^2+(NodeStress(1,I)-NodeStress(3,I))^2+....
        (NodeStress(2,I)-NodeStress(3,I))^2+6*(NodeStress(4,I)^2+NodeStress(5,I)^2+NodeStress(6,I)^2));
end

newNodes=Nodes';
newNodes=newNodes(:);
SF=5.0e2; % ���ηŴ�ϵ��
newNodes=newNodes+SF*U;
newNodes=reshape(newNodes,[3,size(Nodes,1)]);
newNodes=newNodes';
% ����ԭʼ����ͱ����������ͼ
  for i=1:1:size(Elements,1)
    points=Nodes(Elements(i,:),:);
    deformed_points=newNodes(Elements(i,:),:);
    mesh=1:1:4;% ������Ϣ
    % �����嵥Ԫ�ڵ�����
    vertices_matrix = [points(mesh(1,:),1),points(mesh(1,:),2),points(mesh(1,:),3)];
    deformed_vertices_matrix=[deformed_points(mesh(1,:),1),deformed_points(mesh(1,:),2),deformed_points(mesh(1,:),3)];
    % �����嵥Ԫ�ڵ�˳��
    faces_matrix= [1 2 4;2 3 4;3 1 4;1 3 2];% ����ÿ����Ľڵ����
    patch('vertices', vertices_matrix,'faces',faces_matrix,'facecolor','g','FaceAlpha',.5);
    view(3);hold on% ��ͼ
    patch('vertices', deformed_vertices_matrix,'faces',faces_matrix,'facecolor','r','FaceAlpha',.5);
    
  end 
axis equal
view(3);
alpha(1);

% ����Umag��ͼ
Umag=zeros(NodeNum,1);
for i=1:NodeNum
    Umag(i)=sqrt(U(3*i-2)^2+U(3*i-1)^2+U(3*i)^2);
end
PostContour(Nodes,Elements,U,Umag)
title('λ��')
% ����MISES��ͼ
PostContour(Nodes,Elements,U,MISES)
title('Ӧ��')