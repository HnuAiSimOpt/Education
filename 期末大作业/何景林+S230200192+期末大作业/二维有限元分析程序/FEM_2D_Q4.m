%% 注释
%{                 
                                                    有限元期末大作业
  姓名：何景林      学号：S230200192
  1.有限元程序：解决平面应力问题
  2.单元类型：4节点四边形单元
  3.载荷被施加在梁的一侧，而相反的面是固定的
%}
clear
format short
first_time=cputime;             % 用于计算程序执行总时长
%% 初始化
% 物理参数
E=2.1e5;                        %杨氏模量
miu=0.3;                        %泊松比
% 网格节点信息输入
node=importdata('nodes1.txt');  % 节点坐标信息，表示节点x和y方向的坐标
ele=importdata('element1.txt'); % 单元节点信息，表示每个单元上的节点号码
nel=length(ele(:,1));           % 单元总个数
nnode=length(node(:,1));        % 节点总个数
nnel=4;                         % 每单元节点数
ndof=2;                         % 每节点自由度数
sdof=nnode*ndof;                % 总自由度数
edof=nnel*ndof;                 % 每单元自由度数
% 显示网格划分情况
figure(1)
hold on
axis off
axis equal
xp=[];
yp=[];
for ie=1:nel
    for j=1:nnel+1
        j1=mod(j-1,nnel)+1;
        xp(j)=node(ele(ie,j1),1);
        yp(j)=node(ele(ie,j1),2);

    end
    plot(xp,yp,'b-')
end
%% 施加边界条件和载荷
fix_node=[26,  28, 12, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141,142, 143, 144, 145, 146];
bcdof=[];               % 受约束的自由度
bcval=[];               % 受约束自由度对应的位移约束值
for i=1:length(fix_node)             
    bcdof=[bcdof 2*fix_node(i)-1 2*fix_node(i)];
    bcval=[bcval 0 0];
end
% 施加载荷
ff(54,1)=-1000;                             % 给右上角节点施加一个y方向的向下的集中力，大小为1000N

%% 单元刚度矩阵计算与总体刚度矩阵的组装
% 矩阵初始化
ff=zeros(sdof,1);			% 系统载荷初始化
k=zeros(edof,edof);		    % 单元刚度初始化
kk=sparse(sdof,sdof);		% 整体刚度初始化
disp=zeros(sdof,1);         % 整体位移初始化
index=zeros(edof,1);		% 各单元的8个自由度对应的整体的自由度编号
B=sparse(3,edof);           % 应变矩阵B
D=E/(1-miu*miu)*[1  miu 0;
   miu  1  0;
   0  0  (1-miu)/2];        % 平面应力问题的本构矩阵D

G=2;                        % 高斯积分积分点的个数，G=1时为一点高斯积分
[point,weight]=fun_GaussQuadrature(G);
                            % fun_GaussQuadrature()函数返回高斯积分的积分点坐标和权值，具体含义在子函数中有解释
for iel=1:nel               % 遍历各单元
    for i=1:nnel            % 遍历每单元各节点
        nd(i)=ele(iel,i);             % 向量 nd 用来提取每个单元的节点序号
        xcoord(i)=node(nd(i),1);      % 向量 xcoord 用来提取每个单元的节点序号在物理空间下的x坐标
        ycoord(i)=node(nd(i),2);      % 向量 ycoord 用来提取每个单元的节点序号在物理空间下的y坐标
    end
    
    % 开始计算刚度矩阵
    for ix=1:2                      % 遍历x方向的各积分点
        x=point(ix,1);              % 提取x方向积分点的坐标值
        wx=weight(ix,1);            % 提取x方向积分点的权重
        for iy=1:2                  % 遍历y方向的各积分点
            y=point(iy,2);          % 提取y方向积分点的坐标值
            wy=weight(iy,2);        % 提取y方向积分点的权重
            
            [shape,dNdr,dNds]=fun_shapeFunction(x,y);             % 计算积分点的形函数的值和对自然坐标r、s的导数
            
            Jacobian=fun_Jacobian(nnel,dNdr,dNds,xcoord,ycoord);  % 计算雅可比矩阵
            
            detJacobian=det(Jacobian);                            % 计算雅可比矩阵的行列式           
            
            invJacobian=inv(Jacobian);                            % 计算雅可比矩阵的逆
            
            [dNdx,dNdy]=fun_dNdx_dNdy(nnel,dNdr,dNds,invJacobian);% 计算形函数对物理空间坐标x和y的导数
                                                                  % dNdx和dNdy是4维的向量，分别表示第i个形函数对x和y的导数
            
            B=fun_B(nnel,dNdx,dNdy);                              % 组装B矩阵
   
            k=k+B'*D*B*wx*wy*detJacobian;                         % 组装单元刚度矩阵
            
        end
    end
    
    index=fun_index(nd,nnel,ndof);      % 提取单元的系统自由度编号
    
   	kk=fun_assembleKK(kk,k,index);      % 组装整体刚度矩阵
end

%% 零位移边界条件处理：对角线元素改1法
[kk,ff]=fun_boundary(kk,ff,bcdof,bcval);

%% 方程求解
[L,U]=lu(kk);
temp=L\ff;
disp=U\temp ;                               % 得到节点位移

solve_time = cputime-first_time             % 输出求解时间

% 变形后的坐标
disp=load('disp.txt');
disp_x=disp(1:2:sdof);
disp_y=disp(2:2:sdof);
dis=[disp_x,disp_y];
node_new=node+dis;
%% 单元应力应变应变计算
const=[1.866 -0.5 0.134 -0.5;         % 常数矩阵，用于应力外推，用高斯点的应力值来插值节点的应力值
      -0.5 1.866 -0.5 0.134;
      0.134 -0.5 1.866 -0.5;
      -0.5 0.134 -0.5 1.866];
 stress=zeros(nel,4,3);
 strain=zeros(nel,4,3);
for iel=1:nel                         % 遍历各单元
    for i=1:nnel                      % 遍历每单元各节点
        nd(i)=ele(iel,i);             % 向量 nd 用来提取每个单元的节点序号
        xcoord(i)=node(nd(i),1);      % 向量 xcoord 用来提取每个单元的节点序号在物理空间下的x坐标
        ycoord(i)=node(nd(i),2);      % 向量 ycoord 用来提取每个单元的节点序号在物理空间下的y坐标
    end
    pk=0;
    stress_Gauss=zeros(3,4);
    strain_Gauss=zeros(3,4);
    for ix=1:2                        % 遍历x方向的各积分点
        x=point(ix,1);                % 提取x方向积分点的坐标值
        wx=weight(ix,1);              % 提取x方向积分点的权重
        for iy=1:2                    % 遍历y方向的各积分点
            y=point(iy,2);            % 提取y方向积分点的坐标值
            wy=weight(iy,2);          % 提取y方向积分点的权重
            
            [shape,dNdr,dNds]=fun_shapeFunction(x,y);             % 计算积分点的形函数的值和对自然坐标r、s的导数
            
            Jacobian=fun_Jacobian(nnel,dNdr,dNds,xcoord,ycoord);  % 计算雅可比矩阵
            
            detJacobian=det(Jacobian);                            % 计算雅可比矩阵的行列式           
            
            invJacobian=inv(Jacobian);                            % 计算雅可比矩阵的逆
            
            [dNdx,dNdy]=fun_dNdx_dNdy(nnel,dNdr,dNds,invJacobian);% 计算形函数对物理空间坐标x和y的导数
                                                                  % dNdx和dNdy是4维的向量，分别表示第i个形函数对x和y的导数
            
            B=fun_B(nnel,dNdx,dNdy);                              % 组装B矩阵
   
            index=fun_index(nd,nnel,ndof);                        % 单元节点自由度索引
            for i=1:edof
                disp_ele(i,1)=disp(index(i));                     % 单元节点位移
            end            
            strain_ele=B*disp_ele;                                % 计算单元应变
            stress_ele=D*strain_ele;                              % 计算单元应力
            pk=pk+1;
            stress_Gauss(:,pk)=stress_ele;
            strain_Gauss(:,pk)=strain_ele;
        end
    end
    for i=1:3
        for j=1:4
            for n=1:4
                stress(iel,j,i)= stress(iel,j,i)+const(j,n)*stress_Gauss(i,n);
                strain(iel,j,i)= strain(iel,j,i)+const(j,n)*strain_Gauss(i,n);
            end
        end
    end
    
end

% 找出每节点相邻的单元
neigh_node = cell(nnode,1);
neigh_node_ind = cell(nnode,1);
indneigh=zeros(1,nnode);
for i=1:nel
    for j=1:4
        indneigh(ele(i,j))=indneigh(ele(i,j))+1;
        neigh_node{ele(i,j)}(indneigh(ele(i,j)))=i;
        neigh_node_ind{ele(i,j)}(indneigh(ele(i,j)))=j;
    end
end

% 计算节点处的应力应变
stress_node=zeros(3,nnode);
strain_node=zeros(3,nnode);
for inode=1:nnode
    numel= indneigh(inode);
    for i=1:numel
        ind_nel= neigh_node{inode}(i);
        ind_nod=neigh_node_ind{inode}(i);
        for j=1:3
            stress_node(j,inode)=stress_node(j,inode)+stress(ind_nel,ind_nod,j);
            strain_node(j,inode)=strain_node(j,inode)+strain(ind_nel,ind_nod,j);
        end
    end
    stress_node(:,inode)=stress_node(:,inode)/numel;       % 节点处的应力为不同相邻单元磨平后应力的平均值
    strain_node(:,inode)=strain_node(:,inode)/numel;
end

%% 绘图
  fun_picture(node,disp,stress_node,strain_node);
%% 输出
% 计算误差检查
fun_check(disp,stress_node,strain_node);    

% 输出所有计算结果到report_all.txt文件
%-------------------------------------------------------------
fid_out=fopen('report_all.txt','w');
fprintf(fid_out,'TITLE="程序计算结果"\n');
fprintf(fid_out,'VARIABLES="x" "y" "u" "v" "sigax"  "sigmay" "sigmaxy" "epsilonx"  "epsilony" "epsilonxy"\n');
fprintf(fid_out,' N= %8d,E=%8d,Element Type = QUADRILATERA\n',nnode,nel);
for i=1:nnode
    fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',node(i,1),node(i,2),disp(2*i-1),disp(2*i),stress_node(1,i),stress_node(2,i),stress_node(3,i),strain_node(1,i),strain_node(2,i),strain_node(3,i));
end
for i=1:nel
    fprintf(fid_out,'%8d%8d%8d%8d\n',ele(i,1),ele(i,2),ele(i,3),ele(i,4));
end
%% 子函数
% 提取单元的整体自由度编号
function index=fun_index(nd,nnel,ndof)
% nd：单元节点号
% nnel：每单元节点数
% ndof：每节点自由度
 edof = nnel*ndof;
   k=0;
   for i=1:nnel
     start = (nd(i)-1)*ndof;
       for j=1:ndof
         k=k+1;
         index(k)=start+j;
       end
   end

end
% 计算形函数对物理坐标x和y的导数
function [dNdx,dNdy]=fun_dNdx_dNdy(nnel,dNdr,dNds,invJacobian) 
for i=1:nnel
 dNdx(i)=invJacobian(1,1)*dNdr(i)+invJacobian(1,2)*dNds(i);
 dNdy(i)=invJacobian(2,1)*dNdr(i)+invJacobian(2,2)*dNds(i);
end
end
% 返回高斯积分点处的形函数值和形函数对自然坐标的导数值
function [shape,dNdr,dNds]=fun_shapeFunction(r,s)
% r，s为积分点坐标值
% 计算自然坐标系下每个节点的形函数
 shape(1)=0.25*(1-r)*(1-s);
 shape(2)=0.25*(1+r)*(1-s);
 shape(3)=0.25*(1+r)*(1+s);
 shape(4)=0.25*(1-r)*(1+s);
 
% 计算形函数对自然坐标r的导数
 dNdr(1)=-0.25*(1-s);
 dNdr(2)=0.25*(1-s);
 dNdr(3)=0.25*(1+s);
 dNdr(4)=-0.25*(1+s);
 
% 计算形函数对自然坐标s的导数
 dNds(1)=-0.25*(1-r);
 dNds(2)=-0.25*(1+r);
 dNds(3)=0.25*(1+r);
 dNds(4)=0.25*(1-r);
end
% 计算雅可比矩阵
function Jacobian=fun_Jacobian(nnel,dNdr,dNds,xcoord,ycoord)
% nnel:单元节点数
% dNdr，dNds：形函数对自然坐标的导数
% xcoord,ycoord：物理空间下单元四个节点的坐标
 Jacobian=zeros(2,2);
 for i=1:nnel
     Jacobian(1,1)=Jacobian(1,1)+dNdr(i)*xcoord(i);
     Jacobian(1,2)=Jacobian(1,2)+dNdr(i)*ycoord(i);
     Jacobian(2,1)=Jacobian(2,1)+dNds(i)*xcoord(i);
     Jacobian(2,2)=Jacobian(2,2)+dNds(i)*ycoord(i);
 end
end
% 输出二维高斯积分的求积节点
function [point,weight]=fun_GaussQuadrature(i)
% point 和 weight 矩阵的第一列表示x轴积分的节点和权重 
% point 和 weight 矩阵的第二列表示y轴积分的节点和权重
% 有几个求积节点 point 和 weight 就有几行
   point=zeros(i,2);        % i表示有i个求积节点，2指代x轴和y轴
   weight=zeros(i,2);
   
   if i==1                  % 1点高斯积分
       point=[0,0];
       weight=[2,2];
   elseif i==2              % 2点高斯积分
       point=[-0.577,-0.577;
               0.577,0.577];
       weight=[1,1
               1,1];
   end
    
   
end
% 组装整体刚度矩阵
function kk=fun_assembleKK(kk,k,index)
edof = length(index);
for i=1:edof
    ii=index(i);
    for j=1:edof
        jj=index(j);
        kk(ii,jj)=kk(ii,jj)+k(i,j);
    end
end
end
% 组装应变矩阵B
function B=fun_B(nnel,dNdx,dNdy) 
for i=1:nnel
    i1=(i-1)*2+1;  
    i2=i1+1;
    B(1,i1)=dNdx(i);
    B(2,i2)=dNdy(i);
    B(3,i1)=dNdy(i);
    B(3,i2)=dNdx(i);
end
end
% 引入边界条件，采用对角线元素改1法消去矩阵线性方程组的约束项
function [kk,ff]=fun_boundary(kk,ff,bcdof,bcval)
 n=length(bcdof);
 sdof=length(kk);

 for i=1:n
    c=bcdof(i);
    for j=1:sdof
       kk(c,j)=0;
    end

    kk(c,c)=1;
    ff(c)=bcval(i);
 end
end
% 绘制云图
function fun_picture(node,disp,stress_node,strain_node)
x=node(:,1);
y=node(:,2);
% 绘制位移云图
[m,n]=size(disp);
z=[];
for i=1:2:m-1
    z=[z;sqrt(disp(i)*disp(i)+disp(i+1)*disp(i+1))];
end
[X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');  %插值
%去掉模型外部分的值
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z(j,i)=NaN;
        end
    end
end
%绘制等高线图
figure('name','合位移云图');
contourf(X,Y,Z,'LineStyle','none');
colormap jet
hold on;
colorbar;
axis equal
title('合位移云图');

%绘制应力云图
[m1,n2]=size(stress_node);
z2=[];
for i=1:n2
        z2=[z2 sqrt(stress_node(1,i)*stress_node(1,i)+stress_node(2,i)*stress_node(2,i))];
end
[X,Y,Z0]=griddata(x,y,z2,linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear'); 
%去掉模型外部分的值
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z0(j,i)=NaN;
        end
    end
end
%绘制等高线图
figure('name','总应力云图');
contourf(X,Y,Z0,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('总应力云图');

%%x轴应力
[X,Y,Z1]=griddata(x,y,stress_node(1,:),linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');
%去掉模型外部分的值
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z1(j,i)=NaN;
        end
    end
end
%绘制等高线图
figure('name','x轴方向应力云图');
contourf(X,Y,Z1,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('x轴方向应力云图');

%%y轴应力
%插值
[X,Y,Z2]=griddata(x,y,stress_node(2,:),linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');
%去掉模型外部分的值
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z2(j,i)=NaN;
        end
    end
end
%绘制等高线图
figure('name','y轴方向应力云图');
contourf(X,Y,Z2,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('y轴方向应力云图');

%%切应力
[X,Y,Z3]=griddata(x,y,stress_node(3,:),linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');
%去掉模型外部分的值
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z3(j,i)=NaN;
        end
    end
end
%绘制等高线图
figure('name','切应力云图');
contourf(X,Y,Z3,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('切应力云图');

%%绘制应变云图
%%x轴应变
[X,Y,Z1]=griddata(x,y,strain_node(1,:),linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');
%去掉模型外部分的值
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z1(j,i)=NaN;
        end
    end
end
%绘制等高线图
figure('name','x轴方向应变云图');
contourf(X,Y,Z1,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('x轴方向应变云图');

%%y轴应变
%插值
[X,Y,Z2]=griddata(x,y,strain_node(2,:),linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');
%去掉模型外部分的值
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z2(j,i)=NaN;
        end
    end
end
%绘制等高线图
figure('name','y轴方向应变云图');
contourf(X,Y,Z2,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('y轴方向应变云图');

%%切应变
[X,Y,Z3]=griddata(x,y,strain_node(3,:),linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');
%去掉模型外部分的值
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z3(j,i)=NaN;
        end
    end
end
%绘制等高线图
figure('name','切应变云图');
contourf(X,Y,Z3,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('切应变云图');
end
% 与软件对比求解误差
function fun_check(disp,stress,strain)
stress_node=stress';
strain_node=strain';
% 导入位移、应力、应变矩阵
u0=load('u_abaqus.txt');      
[m,n]=size(u0);
for i=1:m
    u(2*i-1, 1)=u0(i, 1);
    u(2*i,1)=u0(i,2);
end
s= load('s_abaqus.txt');
E= load('E_abaqus.txt');
% 计算误差 
d_u=disp-u;         %位移差值
d_s=stress_node-s;  %应力差值
d_E=strain_node-E;  %应变差值

for i=1:2*m
    eu(i,1)=abs(d_u(i,1)/disp(i,1));
end
for j=1:m
    esx(j,1)=abs(d_s(j,1)/stress_node(j,1));
    esy(j,1)=abs(d_s(j,2)/stress_node(j,2));
    esxy(j,1)=abs(d_s(j,3)/stress_node(j,3));
    
    eEx(j,1)=abs(d_E(j,1)/strain_node(j,1));
    eEy(j,1)=abs(d_E(j,2)/strain_node(j,2));
    eExy(j,1)=abs(d_E(j,3)/strain_node(j,3));
end
es=[esx;esy;esxy];
eE=[eEx;eEy;eExy];
% 计算误差均值
eu_total=mean(eu);
es_total=mean(es);
eE_total=mean(eE);
end