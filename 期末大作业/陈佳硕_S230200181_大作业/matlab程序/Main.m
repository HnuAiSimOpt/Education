%%%主程序
%%%用有限元方法进行T形悬臂梁静力学分析（一阶六面体单元）
%%%姓名：陈佳硕
%%%学号：S230200181
clc;clear;
%读取文件获得节点坐标信息Nodes及单元信息Elements
[Nodes, Elements] = Readmessage( 'T.inp' );
Nodes=Nodes*0.001;%单位换算（mm--m）
Mg=5000;%定义T形梁所受载荷（N）
% 外力矩阵 Forces=[受力节点  受力方向(1,2,3分别代表x,y,z)  外力大小]
row_number=6;%受力区域的节点行数
arrange_number=12;%受力区域的节点列数
gap=77;%每行间隔距离
Forces=zeros(row_number*arrange_number,3);
Forces(:,2)=3;%受力方向为z轴方向
Forces(:,3)=-(Mg/(row_number*arrange_number));%每个受力节点的载荷，负号表示z轴负方向
for i=1:row_number
    Forces(1+(i-1)*arrange_number,1)=1+(i-1)*gap;
    Forces(2+(i-1)*arrange_number,1)=14+(i-1)*gap;
    Forces(3+(i-1)*arrange_number:3+(i-1)*arrange_number+9,1)=41+(i-1)*gap:41+(i-1)*gap+9;
end
%约束节点的编号
ConNumber=3851:1:3927;
%约束矩阵 Constraints=[强制位移节点  强制位移方向(1,2,3分别代表x,y,z)  强制位移大小]  
Constraints=zeros(size(ConNumber,2)*3,3);
for i=1:size(ConNumber,2)
Constraints(3*i-2:3*i,:)=[ConNumber(i) 1 0;ConNumber(i) 2 0;ConNumber(i) 3 0;];
end
E=2.1e11; %弹性模量
u=0.3;    %泊松比
%调用函数求解应变应力矩阵D
D=LinearIsotropicD(E,u);
%求解位移U
Dof=3;
NodeCount = size(Nodes,1);  % 节点个数
ElementCount= size(Elements,1);  %单元个数
Dofs = Dof*NodeCount; %总自由度数
U=sparse(Dofs,1);	% 初始化结构位移
K = sparse(Dofs,Dofs);  %初始化总体刚度阵
Force = sparse(Dofs,1);  %初始化外力向量
%计算应力-应变矩阵D
D=LinearIsotropicD(E,u);
for I=1:ElementCount
    % 单元节点坐标
    ElementNodeCoordinate=Nodes(Elements(I,:),:);
    % 计算单刚
    ElementStiffnessMatrix=Ke(D,ElementNodeCoordinate);
    % 计算单元节点自由度编号
    ElementNodeDOF=zeros(1,24);
    for J=1:8
        II=(J-1)*Dof+1;
        ElementNodeDOF(II:II+2)=(Elements(I,J)-1)*Dof+1:(Elements(I,J)-1)*Dof+3;
    end
    K(ElementNodeDOF,ElementNodeDOF)=K(ElementNodeDOF,ElementNodeDOF)+ElementStiffnessMatrix;%根据单元节点自由度  组装总刚
end
% 施加外力
if size(Forces,1)>0
    ForceDOF = Dof*(Forces(:,1)-1)+Forces(:,2);   %计算外力自由度编号
    Force(ForceDOF) = Force(ForceDOF) + Forces(:,3);
end
% 乘大数法施加位移约束
BigNumber=1e8;
ConstraintsNumber=size(Constraints,1);
if ConstraintsNumber~=0
    FixedDof=Dof*(Constraints(:,1)-1)+Constraints(:,2);  %被约束的自由度编号(列向量)
    for i=1:ConstraintsNumber
        K(FixedDof(i),FixedDof(i))=K(FixedDof(i),FixedDof(i))*BigNumber;
        Force(FixedDof(i))=Constraints(i,3)*K(FixedDof(i),FixedDof(i));
    end
end
%计算位移
U = K\Force;
% 求解应力和应变并输出相应云图
OutputResults(Nodes,Elements,D,U)%调用输出结果文件