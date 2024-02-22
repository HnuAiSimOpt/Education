
function [NodeStrain,NodeStress]=NodeStressStrain(U,D,Nodes,Elements)  
ElementNum= size(Elements,1);  %单元个数
Dof=3; 
NodeStrain=zeros(6,ElementNum*4);%11 22 33 12 13 23
NodeStress=zeros(6,ElementNum*4);
Strain = zeros(6,1);
Stress = zeros(6,1);
%循环组装总刚
for i=1:ElementNum
    % 单元节点坐标
    ElementNodeCoordinate=Nodes(Elements(i,:),:); %4*3
    % 计算单元节点自由度编号,即索引位置
    ElementNodeDOF=zeros(1,12);%3*4
    for j=1:4
        ii=(j-1)*Dof+1;
        %单元各节点在总刚度矩阵中的位置索引
        ElementNodeDOF(ii:ii+2)=(Elements(i,j)-1)*Dof+1:(Elements(i,j)-1)*Dof+3;
    end
    % 计算形函数导数
    [NDxyz, ~] = ShapeFunction(ElementNodeCoordinate);%[DN1Dx DN2Dx DN3Dx;DN1Dy DN2Dy DN3Dy;……]
    ElementNodeDisplacement=U(ElementNodeDOF);%12*1 节点位移列阵
    ElementNodeDisplacement=reshape(ElementNodeDisplacement,Dof,4);%3*4
    % 计算单元应变 Strain3_3 3*3的应变矩阵
    Strain3_3=ElementNodeDisplacement*NDxyz';%高斯积分点处应变 3*4  4*3
    %把单元应变矩阵改写成6*1
    Strain=[Strain3_3(1,1) Strain3_3(2,2) Strain3_3(3,3) ...
    Strain3_3(1,2)+Strain3_3(2,1)....
    Strain3_3(2,3)+Strain3_3(3,2) Strain3_3(1,3)+Strain3_3(3,1)]';
    Stress(1:6,1) = D*Strain;%高斯积分点处应变
    %求得节点应力应变
    for X = 1:4
        NodeStrain(1:6,((i-1)*4+X)) = Strain(1:6,1);%按照各单元4节点依次排列1号单元4节点；2号单元4节点
        NodeStress(1:6,((i-1)*4+X)) = Stress(1:6,1);
    end
end
end