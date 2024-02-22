clear
clc

%网格单元数量，m行n列
m=50;
n=20;

%二维矩形几何坐标

vertex=[0,0;0,1;5,1;5,0];


x0=[];
for i=1:m+1
    for j=1:n+1
        x0=[x0; (i-1)*5/m  -1*(1-(j-1)/n)];
    end
end


nodes=[];
for i=1:m
    for j=1:n
        nodes=[nodes; (n+1)*(i-1)+j (n+1)*i+j (n+1)*i+j+1 (n+1)*(i-1)+j+1;];
    end
end

% 绘制网格
iplot=1;
if iplot==1
   figure(1)
   hold on
   axis off
   axis equal
   for ie=1:m*n
        for j=1:4+1
            j1=mod(j-1,4)+1;
            xp(j)=x0(nodes(ie,j1),1);
            yp(j)=x0(nodes(ie,j1),2);
        end
        plot(xp,yp,'-')
   end
          
end
%-------------------------------------------------------------------------
%节点编号
for ii=1:m %第ii行单元
    for jj=1:n %第jj列单元
        node((ii-1)*n+jj,1)=(ii-1)*(n+1)+jj;
        node((ii-1)*n+jj,2)=(ii-1)*(n+1)+jj+1;
        node((ii-1)*n+jj,3)=ii*(n+1)+jj+1;
        node((ii-1)*n+jj,4)=ii*(n+1)+jj;
    end
end
%-------------------------------------------------------------------------
for nd=0:n
    x(nd+1)=vertex(1,1)+nd*(vertex(2,1)-vertex(1,1))/n;   %下边网格节点划分，x坐标
    x((n+1)*m+nd+1)=vertex(4,1)+nd*(vertex(3,1)-vertex(4,1))/n;   %上边网格节点划分，x坐标
    y(nd+1)=vertex(1,2)+nd*(vertex(2,2)-vertex(1,2))/n;   %下边网格节点划分，y坐标
    y((n+1)*m+nd+1)=vertex(4,2)+nd*(vertex(3,2)-vertex(4,2))/n;   %上边网格节点划分，y坐标
end

%内部节点赋值
for ndin=1:n+1  %内部第ndin列节点
    for mdin=1:m-1  %内部第mdin行节点
        x((n+1)*mdin+ndin)=x(ndin)+(x((n+1)*m+ndin)-x(ndin))/m*mdin;   %内部节点横坐标
        y((n+1)*mdin+ndin)=y(ndin)+(y((n+1)*m+ndin)-y(ndin))/m*mdin;   %内部节点纵坐标
    end
end
%--------------------------------------------------------

w1=1;w2=1;%高斯积分权重系数
E=1;%杨氏模量
v=0.3;%泊松比
D=E/(1-v^2)*[1 v 0;v 1 0; 0 0 (1-v)/2];%应变-应力关系转换式

%四个积分点的等参坐标
d=[-1/sqrt(3),-1/sqrt(3);
-1/sqrt(3),1/sqrt(3);
1/sqrt(3),-1/sqrt(3);
1/sqrt(3),1/sqrt(3)];

%四边形各节点Jacob矩阵、应变矩阵、刚度矩阵

%形函数对eta求偏导，eta为自然坐标横坐标
N_eta=zeros(4);
syms eta
for kk=1:4
    zeta=d(kk,2);
    %形函数
    N1=1/4*(1-eta)*(1-zeta);
    N2=1/4*(1+eta)*(1-zeta);
    N3=1/4*(1+eta)*(1+zeta);
    N4=1/4*(1-eta)*(1+zeta);
    N_eta(kk,1:4)=[diff(N1,eta),diff(N2,eta),diff(N3,eta),diff(N4,eta)];
end
%形函数对zeta求偏导，eta为自然坐标纵坐标
N_zeta=zeros(4);
syms zeta
for kk=1:4
    eta=d(kk,1);
    %形函数
    N1=1/4*(1-eta)*(1-zeta);
    N2=1/4*(1+eta)*(1-zeta);
    N3=1/4*(1+eta)*(1+zeta);
    N4=1/4*(1-eta)*(1+zeta);
    N_zeta(kk,1:4)=[diff(N1,zeta),diff(N2,zeta),diff(N3,zeta),diff(N4,zeta)];
end

Ke=zeros((m+1)*(n+1)*2);   %刚度矩阵初始化

%逐个计算各单元刚度矩阵并叠加
for unit=1:m*n
   
    %计算单个单元刚度矩阵
    unitx_loc=[x(node(unit,1)),x(node(unit,2)),x(node(unit,3)),x(node(unit,4))]; 
    unity_loc=[y(node(unit,1)),y(node(unit,2)),y(node(unit,3)),y(node(unit,4))]; 
    unit_loc=[unitx_loc',unity_loc'];   %第unit个单元的物理（几何）坐标
    K=zeros(8,8);   %单个单元刚度矩阵初始化
   for mm=1:4   %mm为节点号（局部，1~4）
      J=[N_eta(mm,:);N_zeta(mm,:)]*unit_loc;   %Jacob矩阵
      Nxy=inv(J)*[N_eta(mm,:);N_zeta(mm,:)];
      B=[Nxy(1,1),0,Nxy(1,2),0,Nxy(1,3),0,Nxy(1,4),0;
      0,Nxy(2,1),0,Nxy(2,2),0,Nxy(2,3),0,Nxy(2,4);
      Nxy(2,1),Nxy(1,1),Nxy(2,2),Nxy(1,2),Nxy(2,3),Nxy(1,3),Nxy(2,4),Nxy(1,4)];   %应变矩阵
      K=K+w1*w2*B'*D*B*det(J);   %刚度矩阵
   end
   
   %单元刚度矩阵叠加
   for p=1:4
       for q=1:4
          Ke(2*node(unit,p)-1:2*node(unit,p),2*node(unit,q)-1:2*node(unit,q))=...
          Ke(2*node(unit,p)-1:2*node(unit,p),2*node(unit,q)-1:2*node(unit,q))+K(2*p-1:2*p,2*q-1:2*q);
       end
   end
end

%施加力矩阵；
Fe=zeros(2*(m+1)*(n+1),1);
Fe(length(Fe),1)=1;

bcdof=[];
bcval=[];
for i=1:21
        bcdof=[bcdof 1+2*(i-1) 2+2*(i-1)];
        bcval=[bcval 0 0];
end

%-----------------------------------------------------------------------
%施加边界条件
 nn=length(bcdof);
 sdof=size(Ke);

 for i=1:nn
    c=bcdof(i);
    for j=1:sdof
       Ke(c,j)=0;
    end

    Ke(c,c)=1;
   Fe(c)=bcval(i);
 end

%-----------------------------------------------------------------------
%计算位移
[LL UU]=lu(Ke);

utemp=LL\Fe;

uv=UU\utemp;



%-----------------------------------------------------------------------
%计算应力矩阵

stress=zeros((m+1)*(n+1),3);   %应力矩阵初始化
a=1+sqrt(3)/2;b=-1/2;c=1-sqrt(3)/2;   %用于局部应力矩阵到整体应力矩阵转换
for unit=1:m*n

    unitx_loc=[x(node(unit,1)),x(node(unit,2)),x(node(unit,3)),x(node(unit,4))]; 
    unity_loc=[y(node(unit,1)),y(node(unit,2)),y(node(unit,3)),y(node(unit,4))]; 
    unit_loc=[unitx_loc',unity_loc'];   %第unit个单元的物理（几何）坐标
    %计算单个单元局部应力矩阵
   for mm=1:4   %mm为节点号（局部，1~4）
      J=[N_eta(mm,:);N_zeta(mm,:)]*unit_loc;   %Jacob矩阵
      Nxy=inv(J)*[N_eta(mm,:);N_zeta(mm,:)];
      B=[Nxy(1,1),0,Nxy(1,2),0,Nxy(1,3),0,Nxy(1,4),0;
      0,Nxy(2,1),0,Nxy(2,2),0,Nxy(2,3),0,Nxy(2,4);
      Nxy(2,1),Nxy(1,1),Nxy(2,2),Nxy(1,2),Nxy(2,3),Nxy(1,3),Nxy(2,4),Nxy(1,4)];   %应变矩阵
      stress_e(mm,:)=[D*B*[uv(2*node(unit,1)-1),uv(2*node(unit,1)),...
      uv(2*node(unit,2)-1),uv(2*node(unit,2)),uv(2*node(unit,3)-1),...
      uv(2*node(unit,3)),uv(2*node(unit,4)-1),uv(2*node(unit,4))]']';   %单元局部应力矩阵
      
   end
   %局部应力矩阵转换整体应力矩阵，并叠加
   stress_whole=[a b c b;b a b c;c b a b;b c b a]*stress_e;
   for mm=1:4
       stress(node(unit,mm),:)=stress(node(unit,mm),:)+stress_whole(mm,:);   %总体应力矩阵叠加
   end
end
stress;
%多单元公用节点应力取平均值
for nod=1:(m+1)*(n+1)
    for row=0:m
    if nod>(n+1)*row+1 & nod<(n+1)*(row+1)
        stress(nod,:)=stress(nod,:)/2;
    end
    end
    for column=0:n
    if nod>column+1 & nod<(n+1)*m+column+1 & (nod-column-1)/(n+1)==fix((nod-column-1)/(n+1))
        stress(nod,:)=stress(nod,:)/2;
    end
    end
end
uv;
stress;
dispu=zeros(1,1071);%x方向位移
dispv=zeros(1,1071);%y方向位移

for i=1:1071
    dispu(i)=uv(2*i-1);
    dispv(i)=uv(2*i);
end


[X,Y,Z]=griddata(x0(:,1),x0(:,2),dispu,linspace(0,5,50)',linspace(0,-1,20),'v4');
figure(2);
pcolor(X,Y,Z);
shading interp; 
colormap(jet); %x方向位移云图
[X,Y,Z]=griddata(x0(:,1),x0(:,2),dispv,linspace(0,5,50)',linspace(0,-1,20),'v4');
figure(3);
pcolor(X,Y,Z);
shading interp; 
colormap(jet); %y方向位移云图
[X,Y,Z]=griddata(x0(:,1),x0(:,2),stress(:,1),linspace(0,5,50)',linspace(0,-1,20),'v4');
figure(4);
pcolor(X,Y,Z);
shading interp; 
colormap(jet); %x方向应力云图
[X,Y,Z]=griddata(x0(:,1),x0(:,2),stress(:,2),linspace(0,5,50)',linspace(0,-1,20),'v4');
figure(5);
pcolor(X,Y,Z);
shading interp; 
colormap(jet); %y方向应力云图
figure(6);
quiver(x0(:,1),x0(:,2),dispu',dispv');%变形图

% 数据输出
%-------------------------------------------------------------
fid_out=fopen('result.txt','w');

fprintf(fid_out,'VARIABLES="x" "y" "u" "v" "sigax"  "sigmay" "sigmaxy"\n');


for i=1:nod
    fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',x0(i,1),x0(i,2),dispu(i)',dispv(i)',stress(i,1),stress(i,2),stress(i,3));
end







