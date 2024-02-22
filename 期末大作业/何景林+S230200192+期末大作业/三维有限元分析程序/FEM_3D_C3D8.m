%% 注释
%{                 
                                                    有限元期末大作业
  姓名：何景林      学号：S230200192
  1.有限元程序：三维悬臂梁的有限元分析
  2.单元类型：C3D8单元
  3.载荷被施加在梁的一侧，而相反的面是固定的
%}
clc
clear
close all
%% 初始化
% 读取节点坐标
Nodes = xlsread('Nodes1.xlsx');
[N,~] = size(Nodes);  %节点数
%  读取单元与节点之间的对应关系   
Elems = xlsread('Elements1.xlsx');
[E,~] = size(Elems);  %单元数
%  每个单元的节点数目
NE = 8;  
%  读取材料信息
Mats = load('Materials.txt');%材料属性
ipstrn = 2;
nstrn = 3;
%% 定义力和边界节点
%  边界约束节点
j_dbc=1;
for (i=1:N)
    if (Nodes(i,4)==0)
        DBC(j_dbc,1)=Nodes(i,1);
        DBC(j_dbc,2)=1;
        DBC(j_dbc,3)=0;
        j_dbc=j_dbc+1;
    end
end
[P,~] = size(DBC);
%  载荷施加节点
j_nbc=1;
for (i=1:N)
    if (Nodes(i,4)==60 && Nodes(i,3)==20)
        right(j_nbc,1)=Nodes(i,1);
        j_nbc=j_nbc+1;
    end
end
j_nbc=1;
for i=1:E
    for j=1:size(right(:,1))
        for k=3:10
            if Elems(i,k)==right(j,1)
                el_list(j_nbc,1)=Elems(i,1);
                el_list(j_nbc,2)=right(j,1);
                j_nbc=j_nbc+1;
                break
            end
        end
    end
end
NBC(:,1)=unique(el_list(:,1));
j_nbc=1;
for i=1:2:size(el_list(:,1))
    for j=3:10
        if (el_list(i,2)==Elems(el_list(i,1),j))
            
            NBC(j_nbc,2)=el_list(i,2);
            NBC(j_nbc,3)=el_list(i+1,2);
            j_nbc=j_nbc+1;
        end
    end
end
NBC(:,4)=2;
NBC(:,5)=1;
[Q,l] = size(NBC); 
%% 计算刚度矩阵
%  总自由度
udof = 3;     % 单个节点自由度
NDOF = N*udof;
%  矩阵初始化
K = zeros(NDOF,NDOF);   % 刚度矩阵
U = zeros(NDOF,1);      % 位移矩阵
F = zeros(NDOF,1);      % 力矩阵
%  位移约束惩罚
Klarge = 10^8;
%  设置高斯点位置和权重
NG = 8; % 高斯点个数
[XG,WG] = Gauss_Points(NG);
%  遍历所有单元
for e = 1:E    
%   建立单元和节点坐标的联系
    Nnums = Elems(e,3:2+NE);
    xyz = Nodes(Nnums(:),2:4);  
%   提取单元弹性杨氏模量和泊松比
    Y = Mats(Elems(e,2),2);
    nu = Mats(Elems(e,2),3);    
%   构建单元刚度矩阵
    [Ke] = Stiff(ipstrn,xyz,Y,nu,udof,NE,NG,XG,WG);
    % 组装整体刚度矩阵
    ig = udof*(Nnums(:)-1);
    for ni = 1:NE
        i0 = udof*(ni-1);
        for nj = 1:NE
            j0 = udof*(nj-1);
            for i = 1:udof
                for j = 1:udof
                    K(ig(ni)+i,ig(nj)+j) = K(ig(ni)+i,ig(nj)+j) + Ke(i0+i,j0+j);
                end
            end
        end
    end
end
%% 定义力矢
NES = 2; %受力单元节点个数
NGS = 2; %受力单元高斯点个数
[XGS,WGS] = Gauss_Points_Surf(NGS);%加载节点高斯点坐标和权重
for q = 1:Q   
    in   = zeros(NES);
    tval = zeros(NES,1);
    fval = zeros(NES,1);
    % 确定加载的节点
    e = NBC(q,1);%加载单元
    in1 = NBC(q,2);%加载单元节点1
    in2 = NBC(q,3);%加载单元节点2
    idof = NBC(q,4);%加载方向
    tval(:,1) = NBC(q,5);  
    for i=1:NGS      
        xi  = XGS(i);
        wgt = WGS(i);
        NshapeS(1) = (1-xi)/2;
        NshapeS(2) = (1+xi)/2;
        DNshapeS(1) = -1/2;
        DNshapeS(2) = +1/2;
        xyS(1,1) = Nodes(in1,2);
        xyS(1,2) = Nodes(in1,3);
        xyS(1,3) = Nodes(in1,4);
        xyS(2,1) = Nodes(in2,2);
        xyS(2,2) = Nodes(in2,3);
        xyS(2,3) = Nodes(in2,4);
        [detJS] = Jacobian_Surf(NES,xi,xyS,DNshapeS); 
        fval = fval + wgt*NshapeS'*NshapeS*tval*detJS;     
    end  
    iloc1 = udof*(in1-1)+idof;
    iloc2 = udof*(in2-1)+idof;     
    F(iloc1) = F(iloc1) + fval(1);
    F(iloc2) = F(iloc2) + fval(2);       
end
%% 定义边界条件
for p = 1:P
    inode = DBC(p,1);
    idof1 =1;
    idof2 =2;
    idof3 =3;
    idiag1 = udof*(inode-1) + idof1;
    idiag2 = udof*(inode-1) + idof2;
    idiag3 = udof*(inode-1) + idof3;
    K(idiag1,idiag1) = Klarge;
    K(idiag2,idiag2) = Klarge;
    K(idiag3,idiag3) = Klarge;
    F(idiag1) = Klarge*DBC(p,3);
    F(idiag2) = Klarge*DBC(p,3);
    F(idiag3) = Klarge*DBC(p,3);
end
 F=F/sum(F);
 %
 %F = F*5;
 F(956)=0.1;
 F(959)=0.1;
 F(995)=0.2;
 F(998)=0.2;
 F(1001)=0.2;
 F(1004)=0.2;
 %}
%% 求解位移、应变和应力
  U = inv(K)*F;

% 初始化位移、应变、应力
nedof = udof*NE;
Disp = zeros(E,nedof);
Eps = zeros(E,nstrn,NG);
Sig = zeros(E,nstrn,NG);
% 计算应变应力
for e = 1:E
    Nnums = Elems(e,3:2+NE);
    xyz = Nodes(Nnums(:),2:4);
    h = Elems(e,3);    
    % 提取单元弹性杨氏模量和泊松比
    Y = Mats(Elems(e,2),2);
    nu = Mats(Elems(e,2),3);
    % 提取位移
    for i=1:NE
        inode = Nnums(i);
        iglb1 = udof*(inode-1)+1;
        iglb2 = udof*inode;
        iloc1 = udof*(i-1)+1;
        iloc2 = udof*i;
        Disp(e,iloc1) = U(iglb1);
        Disp(e,iloc2) = U(iglb2);
    end
    % 应变、应力
    u = Disp(e,:)';
    [eps,sig] = Str(ipstrn,xyz,u,h,Y,nu,udof,NE,NG,XG);
    Eps(e,:,:) = eps(:,:);
    Sig(e,:,:) = sig(:,:);    
end 
%% 绘图
Plot_mesh(Nodes(:,2:4),Elems(:,3:10));
title('未施加载荷');
xlabel('X ');
ylabel('Y ');
zlabel('Z ');
j=1;
for i=1:3:size(U)
n_disp(j,1)=U(i);
n_disp(j,2)=U(i+1);
n_disp(j,3)=U(i+2);
j=j+1;
end
n_final(:,1)=Nodes(:,2)+n_disp(:,1);
n_final(:,2)=Nodes(:,3)+n_disp(:,2);
n_final(:,3)=Nodes(:,4)+n_disp(:,3);
figure;
Plot_mesh(n_final,Elems(:,3:10));
title('施加载荷');
xlabel('X ');
ylabel('Y ');
zlabel('Z ');
%% 子函数
% 计算高斯积分点等参坐标和权重
function [XG,WG] = Gauss_Points(NG)
alf = sqrt(1/3);
XG(1,1) = -alf;
XG(2,1) = +alf;
XG(3,1) = +alf;
XG(4,1) = -alf;
XG(5,1) = -alf;
XG(6,1) = +alf;
XG(7,1) = +alf;
XG(8,1) = -alf;

XG(1,2) = -alf;
XG(2,2) = -alf;
XG(3,2) = +alf;
XG(4,2) = +alf;
XG(5,2) = -alf;
XG(6,2) = -alf;
XG(7,2) = +alf;
XG(8,2) = +alf;

XG(1,3) = -alf;
XG(2,3) = -alf;
XG(3,3) = -alf;
XG(4,3) = -alf;
XG(5,3) = +alf;
XG(6,3) = +alf;
XG(7,3) = +alf;
XG(8,3) = +alf;
for i=1:NG
    WG(i) = 1;
end

end

% 计算单元刚度矩阵
function [Ke] = Stiff(ipstrn,xyz,Y,nu,udof,NE,NG,XG,WG)
ndof = NE*udof;
nstrn = 6;
Ke = zeros(ndof,ndof);

for i=1:NG
    
    xi  = XG(i,1);
    eta = XG(i,2);
    mu = XG(i,3);
    wgt = WG(i);
    % 计算形函数对局部坐标的偏导
    DNshape(1,1)=-((eta - 1)*(mu - 1))/8;
    DNshape(2,1)= ((eta - 1)*(mu - 1))/8;
    DNshape(3,1)=-((eta + 1)*(mu - 1))/8;
    DNshape(4,1)= ((eta + 1)*(mu - 1))/8;
    DNshape(5,1)= ((eta - 1)*(mu + 1))/8;
    DNshape(6,1)=-((eta - 1)*(mu + 1))/8;
    DNshape(7,1)= ((eta + 1)*(mu + 1))/8;
    DNshape(8,1)=-((eta + 1)*(mu + 1))/8;
    
    DNshape(1,2)=-(xi/8 - 1/8)*(mu - 1);
    DNshape(2,2)= (xi/8 + 1/8)*(mu - 1);
    DNshape(3,2)=-(xi/8 + 1/8)*(mu - 1);
    DNshape(4,2)= (xi/8 - 1/8)*(mu - 1);
    DNshape(5,2)= (xi/8 - 1/8)*(mu + 1);
    DNshape(6,2)=-(xi/8 + 1/8)*(mu + 1);
    DNshape(7,2)= (xi/8 + 1/8)*(mu + 1);
    DNshape(8,2)=-(xi/8 - 1/8)*(mu + 1);
    
    DNshape(1,3)=-(xi/8 - 1/8)*(eta - 1);
    DNshape(2,3)= (xi/8 + 1/8)*(eta - 1);
    DNshape(3,3)=-(xi/8 + 1/8)*(eta + 1);
    DNshape(4,3)= (xi/8 - 1/8)*(eta + 1);
    DNshape(5,3)= (xi/8 - 1/8)*(eta - 1);
    DNshape(6,3)=-(xi/8 + 1/8)*(eta - 1);
    DNshape(7,3)= (xi/8 + 1/8)*(eta + 1);
    DNshape(8,3)=-(xi/8 - 1/8)*(eta + 1);
    % 计算雅可比矩阵
    Jac = zeros(3);
    for ii=1:NE
        Jac(1,1) = Jac(1,1) + DNshape(ii,1)*xyz(ii,1);
        Jac(1,2) = Jac(1,2) + DNshape(ii,1)*xyz(ii,2);
        Jac(1,3) = Jac(1,3) + DNshape(ii,1)*xyz(ii,3);
        Jac(2,1) = Jac(2,1) + DNshape(ii,2)*xyz(ii,1);
        Jac(2,2) = Jac(2,2) + DNshape(ii,2)*xyz(ii,2);
        Jac(2,3) = Jac(2,3) + DNshape(ii,2)*xyz(ii,3);
        Jac(3,1) = Jac(3,1) + DNshape(ii,3)*xyz(ii,1);
        Jac(3,2) = Jac(3,2) + DNshape(ii,3)*xyz(ii,2);
        Jac(3,3) = Jac(3,3) + DNshape(ii,3)*xyz(ii,3);
    end
    
    detJ = det(Jac);
    Jhat = inv(Jac);
    %计算B矩阵
    B = zeros(nstrn,ndof);
    i=1;
    for j=1:NE
        qx=DNshape(j,1)*Jhat(1,1)+DNshape(j,2)*Jhat(1,2)+DNshape(j,3)*Jhat(1,3);
        qy=DNshape(j,1)*Jhat(2,1)+DNshape(j,2)*Jhat(2,2)+DNshape(j,3)*Jhat(2,3);
        qz=DNshape(j,1)*Jhat(3,1)+DNshape(j,2)*Jhat(3,2)+DNshape(j,3)*Jhat(3,3);
        
        B(1,i)=qx;
        B(1,i+1)=0;
        B(1,i+2)=0;
        B(2,i)=0;
        B(2,i+1)=qy;
        B(2,i+2)=0;
        B(3,i)=0;
        B(3,i+1)=0;
        B(3,i+2)=qz;
        B(4,i)=qy;
        B(4,i+1)=qx;
        B(4,i+2)=0;
        B(5,i)=0;
        B(5,i+1)=qz;
        B(5,i+2)=qy;
        B(6,i)=qz;
        B(6,i+1)=0;
        B(6,i+2)=qx;
        i=i+3;
    end
    % 计算三维弹性矩阵C
    lambda=nu*Y/((1+nu)*(1-2*nu));
    c=Y/(2*(1+nu));
    C = [lambda+2*nu lambda lambda 0 0 0; lambda lambda+2*nu lambda 0 0 0; lambda lambda lambda+2*nu 0 0 0; 0 0 0 nu 0 0; 0 0 0 0 nu 0; 0 0 0 0 0 nu];
    % 计算单元刚度矩阵
    Ke = Ke + wgt*transpose(B)*C*B*detJ;
    
end

end
% 计算受力处高斯点坐标和权重
function [XGS,WGS] = Gauss_Points_Surf(NGS)

if (NGS == 2)
    
    alf = sqrt(1/3);

    XGS(1,1) = -alf;
    XGS(2,1) = +alf;

    WGS(1) = 1;
    WGS(2) = 1;
    
elseif (NGS == 3)
    
    alf = sqrt(3/5);

    XGS(1,1) = -alf;
    XGS(2,1) = 0;
    XGS(3,1) = +alf;

    WGS(1) = 5/9;
    WGS(2) = 8/9;
    WGS(3) = 5/9;
elseif (NGS == 4)
    alf = 0.8611363115940526;
    bet = 0.3399810435848563;
    
    XGS(1,1) = -alf;
    XGS(2,1) = -bet;
    XGS(3,1) = bet;
    XGS(4,1) = alf;
    
    WGS(1)=0.3478548451374538;
    WGS(2)=0.6521451548625461;
    WGS(3)=0.6521451548625461;
    WGS(4)=0.3478548451374538;
end
end

function [detJS] = Jacobian_Surf(NES,~,xyS,DNshapeS)

dxdxi = 0;
dydxi = 0;
dzdxi = 0;
for i=1:NES
    dxdxi = dxdxi + DNshapeS(i)*xyS(i,1);
    dydxi = dydxi + DNshapeS(i)*xyS(i,2);
    dzdxi = dzdxi + DNshapeS(i)*xyS(i,3);
end

detJS = sqrt( dxdxi*dxdxi + dydxi*dydxi + dzdxi*dzdxi);
end
% 计算应力应变
function [eps,sig] = Str(ipstrn,xyz,u,h ,Y,nu,udof,NE,NG,XG)
ndof = NE*udof;
nstrn = 3;
eps = zeros(nstrn,NG);
sig = zeros(nstrn,NG);

for i=1:NG
    
   xi  = XG(i,1);
   eta = XG(i,2);
   mu = XG(i,3);
   
% 计算形函数对局部坐标的偏导
    DNshape(1,1)=-((eta - 1)*(mu - 1))/8;
    DNshape(2,1)= ((eta - 1)*(mu - 1))/8;
    DNshape(3,1)=-((eta + 1)*(mu - 1))/8;
    DNshape(4,1)= ((eta + 1)*(mu - 1))/8;
    DNshape(5,1)= ((eta - 1)*(mu + 1))/8;
    DNshape(6,1)=-((eta - 1)*(mu + 1))/8;
    DNshape(7,1)= ((eta + 1)*(mu + 1))/8;
    DNshape(8,1)=-((eta + 1)*(mu + 1))/8;
    
    DNshape(1,2)=-(xi/8 - 1/8)*(mu - 1);
    DNshape(2,2)= (xi/8 + 1/8)*(mu - 1);
    DNshape(3,2)=-(xi/8 + 1/8)*(mu - 1);
    DNshape(4,2)= (xi/8 - 1/8)*(mu - 1);
    DNshape(5,2)= (xi/8 - 1/8)*(mu + 1);
    DNshape(6,2)=-(xi/8 + 1/8)*(mu + 1);
    DNshape(7,2)= (xi/8 + 1/8)*(mu + 1);
    DNshape(8,2)=-(xi/8 - 1/8)*(mu + 1);
    
    DNshape(1,3)=-(xi/8 - 1/8)*(eta - 1);
    DNshape(2,3)= (xi/8 + 1/8)*(eta - 1);
    DNshape(3,3)=-(xi/8 + 1/8)*(eta + 1);
    DNshape(4,3)= (xi/8 - 1/8)*(eta + 1);
    DNshape(5,3)= (xi/8 - 1/8)*(eta - 1);
    DNshape(6,3)=-(xi/8 + 1/8)*(eta - 1);
    DNshape(7,3)= (xi/8 + 1/8)*(eta + 1);
    DNshape(8,3)=-(xi/8 - 1/8)*(eta + 1);
    % 计算雅可比矩阵
    Jac = zeros(3);
    for ii=1:NE
        Jac(1,1) = Jac(1,1) + DNshape(ii,1)*xyz(ii,1);
        Jac(1,2) = Jac(1,2) + DNshape(ii,1)*xyz(ii,2);
        Jac(1,3) = Jac(1,3) + DNshape(ii,1)*xyz(ii,3);
        Jac(2,1) = Jac(2,1) + DNshape(ii,2)*xyz(ii,1);
        Jac(2,2) = Jac(2,2) + DNshape(ii,2)*xyz(ii,2);
        Jac(2,3) = Jac(2,3) + DNshape(ii,2)*xyz(ii,3);
        Jac(3,1) = Jac(3,1) + DNshape(ii,3)*xyz(ii,1);
        Jac(3,2) = Jac(3,2) + DNshape(ii,3)*xyz(ii,2);
        Jac(3,3) = Jac(3,3) + DNshape(ii,3)*xyz(ii,3);
    end
    
    detJ = det(Jac);
    Jhat = inv(Jac);
   
   B = zeros(nstrn,ndof);
   for j=1:NE
       jloc1 = 2*(j-1)+1;
       jloc2 = jloc1 + 1;
       B(1,jloc1) = B(1,jloc1) + Jhat(1,1)*DNshape(j,1) ...
           + Jhat(1,2)*DNshape(j,2);
       B(2,jloc2) = B(2,jloc2) + Jhat(2,1)*DNshape(j,1) ...
           + Jhat(2,2)*DNshape(j,2);
       B(3,jloc1) = B(3,jloc1) + Jhat(2,1)*DNshape(j,1) ...
           + Jhat(2,2)*DNshape(j,2);
       B(3,jloc2) = B(3,jloc2) + Jhat(1,1)*DNshape(j,1) ...
           + Jhat(1,2)*DNshape(j,2);
   end
   
   if (ipstrn == 1)
       c = Y*(1-nu)/(1-2*nu)/(1+nu);
       C = c*[ 1 nu/(1-nu) 0; nu/(1-nu) 1 0; 0 0 (1-2*nu)/(1-nu)/2 ];
   else
       c = Y/(1-nu)/(1+nu);
       C = c*[ 1 nu 0; nu 1 0; 0 0 (1-nu)/2 ];
   end
   
   eps(:,i) = B*u;
   sig(:,i) = C*eps(:,i);
   
end
end
% 绘图函数
function Plot_mesh(node_coord,elements)

n_el = length(elements) ;                  % 单元数
node_face = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8]; 
XYZ = cell(1,n_el) ;

for e=1:n_el
    nd=elements(e,:);
    XYZ{e} = [node_coord(nd,1)  node_coord(nd,2) node_coord(nd,3)] ;
end

% Plot 
axis equal;
axis tight;
cellfun(@patch,repmat({'Vertices'},1,n_el),XYZ,repmat({'Faces'},1,n_el),repmat({node_face},1,n_el),repmat({'FaceColor'},1,n_el),repmat({'w'},1,n_el));
        view(3)
end
