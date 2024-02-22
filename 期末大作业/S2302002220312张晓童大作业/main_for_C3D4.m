
clear all;
clc;
close all;
E=210000; %弹性模量
u=0.3;    %泊松比
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
Dof=3;                          % 每个节点三个自由度
NodeNum = size(Nodes,1);        % 节点个数
ElementNum= size(Elements,1);   % 单元个数
Dofs = Dof*NodeNum;             % 总自由度数
U=sparse(Dofs,1);	            % 初始化结构位移（空矩阵），利用稀疏矩阵提高计算速度
K = sparse(Dofs,Dofs);          % 初始化总体刚度阵（空矩阵）
F = sparse(Dofs,1);             % 初始化外力向量（空矩阵）

for I=1:ElementNum
    % 单元节点坐标
    ElementNodeCoordinate=Nodes(Elements(I,:),:);
    % 计算单元刚度矩阵
    ElementStiffnessMatrix=Ke_C3D4(D,ElementNodeCoordinate);
    % 计算单元节点自由度编号
    ElementNodeDOF=zeros(1,12);
    for J=1:4
        II=(J-1)*Dof+1;
        ElementNodeDOF(II:II+2)=(Elements(I,J)-1)*Dof+1:(Elements(I,J)-1)*Dof+3;
    end
    %整体刚度矩阵组装，根据单元节点自由度组装总刚
    K(ElementNodeDOF,ElementNodeDOF)=K(ElementNodeDOF,ElementNodeDOF)+ElementStiffnessMatrix;
end
% 施加外力
if size(Forces,1)>0
    %计算外力自由度编号,第一列为节点号，第二列为方向编号
    ForceDOF = Dof*(Forces(:,1)-1)+Forces(:,2);
    F(ForceDOF) = F(ForceDOF) + Forces(:,3);
end
% 乘大数法（罚函数法）施加位移约束
BigNumber=1e8;
ConsNum=size(Constraints,1);
if ConsNum~=0
    FixedDof=Dof*(Constraints(:,1)-1)+Constraints(:,2);  %被约束的自由度编号(列向量)
    for i=1:ConsNum
        K(FixedDof(i),FixedDof(i))=K(FixedDof(i),FixedDof(i))*BigNumber;
        F(FixedDof(i))=Constraints(i,3)*K(FixedDof(i),FixedDof(i));
    end
end
% 计算位移
U = K\F;
% 输出结果
% 计算高斯点、节点应力应变值
[NodeStrain,NodeStress]=NodeStressStrain(U,D,Nodes,Elements);
% 求得MISES应力矩阵
ElementNodeCount=4; % 每个单元节点数
MISES=zeros(1,ElementNum*ElementNodeCount);
for I=1:ElementNum*ElementNodeCount
    MISES(I)=sqrt(0.5)*sqrt((NodeStress(1,I)-NodeStress(2,I))^2+(NodeStress(1,I)-NodeStress(3,I))^2+....
        (NodeStress(2,I)-NodeStress(3,I))^2+6*(NodeStress(4,I)^2+NodeStress(5,I)^2+NodeStress(6,I)^2));
end

newNodes=Nodes';
newNodes=newNodes(:);
SF=5.0e2; % 变形放大系数
newNodes=newNodes+SF*U;
newNodes=reshape(newNodes,[3,size(Nodes,1)]);
newNodes=newNodes';
% 绘制原始网格和变形网格叠加图
  for i=1:1:size(Elements,1)
    points=Nodes(Elements(i,:),:);
    deformed_points=newNodes(Elements(i,:),:);
    mesh=1:1:4;% 网格信息
    % 四面体单元节点坐标
    vertices_matrix = [points(mesh(1,:),1),points(mesh(1,:),2),points(mesh(1,:),3)];
    deformed_vertices_matrix=[deformed_points(mesh(1,:),1),deformed_points(mesh(1,:),2),deformed_points(mesh(1,:),3)];
    % 四面体单元节点顺序
    faces_matrix= [1 2 4;2 3 4;3 1 4;1 3 2];% 给出每个面的节点序号
    patch('vertices', vertices_matrix,'faces',faces_matrix,'facecolor','g','FaceAlpha',.5);
    view(3);hold on% 绘图
    patch('vertices', deformed_vertices_matrix,'faces',faces_matrix,'facecolor','r','FaceAlpha',.5);
    
  end 
axis equal
view(3);
alpha(1);

% 绘制Umag云图
Umag=zeros(NodeNum,1);
for i=1:NodeNum
    Umag(i)=sqrt(U(3*i-2)^2+U(3*i-1)^2+U(3*i)^2);
end
PostContour(Nodes,Elements,U,Umag)
title('位移')
% 绘制MISES云图
PostContour(Nodes,Elements,U,MISES)
title('应力')