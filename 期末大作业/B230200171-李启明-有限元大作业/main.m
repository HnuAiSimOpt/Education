%%---------------空间八节点六面体单元有限元分析--------------------
clear
clc
close all

%% 材料、几何属性
E=1.0e10;%弹性模量
nu= 0.25;%泊松比

Lx=1;
Ly = 0.1;
Lz = 0.1;
nx = 40;
ny = 4;
nz = 4;


%% 生成八节点六面体网格
[Node,Element] =generate_mesh(Lx,Ly,Lz,nx,ny,nz);
femesh3D(Element,Node);

%节点、单元、自由度个数
nnode = size(Node,1);
nele = size(Element,1);
nDOF = nnode*3;

%% 组装单元刚度矩阵
KK= zeros(nDOF);
for iele=1:nele
    element_iele = Element(iele,:);
    coord = Node(element_iele,:);
    %单元刚度矩阵
    Ke = stiffnessMatrix(coord,E,nu);

    %组装
    DOF = zeros(1,24);
    for i = 1:length(element_iele)
        DOF((i-1)*3+1:(i-1)*3+3) = element_iele(i)*3-2:element_iele(i)*3; 
    end

    for i = 1:length(DOF)
        for j = 1:length(DOF)
            KK(DOF(i),DOF(j)) =  KK(DOF(i),DOF(j)) + Ke(i,j);
        end
    end

end

%% 生成节点载荷列阵
F= zeros(nDOF,1);%初始化载荷列阵
%加载点编号
loadPoint1 = find(abs(Node(:,1)-Lx)<1e-6&abs(Node(:,2)-0)<1e-6&abs(Node(:,3)-Lz)<1e-6);
F(loadPoint1*3) = -1.0e5; %加1e5N向下的力

loadPoint2 = find(abs(Node(:,1)-Lx)<1e-6&abs(Node(:,2)-Ly)<1e-6&abs(Node(:,3)-Lz)<1e-6);
F(loadPoint2*3) = -1.0e5; %加1e5N向下的力

%% 施加边界约束
%被约束的节点编号
node_BC = find(abs(Node(:,1))<1e-6);
BC = [];%被约束的自由度编号
for i = 1:length(node_BC)
        BC((i-1)*3+1:(i-1)*3+3) = node_BC(i)*3-2:node_BC(i)*3; 
end

fdof=[1:nDOF]';%将所有自由度编号从1到最大排成一个矢量
fdof(BC)=[];%将被约束的自由度全部删掉
KKc=KK(fdof,fdof);%保留刚度矩阵中未被约束的行/列
Fc=F(fdof,:);%保留节点载荷列阵中未被约束的元素

%% 求解
uc=KKc\Fc;
%还原为所有自由度的位移（被约束的自由度补0）
U= zeros(nDOF,1);
U(fdof,:)=uc; 
%X方向位移
fprintf('X方向位移：\n')
Ux = U(1:3:end-2)
%Y方向位移
fprintf('Y方向位移：\n')
Uy = U(2:3:end-1)
%Z方向位移
fprintf('Z方向位移：\n')
Uz = U(3:3:end)

%% 计算单元中心的应力，每个单元六个应力分量，stress的行数为单元个数，列数为6
stress = zeros(nele,6);

for iele = 1:nele
    element_iele = Element(iele,:);
    coord = Node(element_iele,:);
     DOF = zeros(1,24);
    for i = 1:length(element_iele)
        DOF((i-1)*3+1:(i-1)*3+3) = element_iele(i)*3-2:element_iele(i)*3; 
    end
    %单元位移Ue
    Ue = U(DOF);
    
    %几何矩阵
    B = MatrixB(coord, 0,0,0);
    %弹性矩阵
    D = MatrixD(E,nu);
    
    stress(iele,:) = D*B*Ue;
    
end


%% 后处理
%绘制x方向位移云图
figure
fecontour3D(Element,Node,Ux);
title('位移Ux','FontSize',18)

%绘制y方向位移云图
figure
fecontour3D(Element,Node,Uy);
title('位移Uy','FontSize',18)

%绘制z方向位移云图
figure
fecontour3D(Element,Node,Uz);
title('位移Uz','FontSize',18)

%应力分量排列：[sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_yz,sigma_xz]
str_stress = {'xx','yy','zz','xy','yz','xz'};

for i = 1:6
    sigma = stress(:,i);
    [Node1,Element1,sigma1]= smooth_stress(sigma,Node,Element);

    figure
    fecontour3D(Element1,Node1,sigma1);
    title(['应力\sigma_{',str_stress{i},'}'],'FontSize',18)

end







