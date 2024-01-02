%%%%%%%%%
%有限元编程作业
%作者：于东洋  
%学号：S230200187
%完成时间：2023/12/25
%程序目的：对一个长X宽为2mX1m的悬臂梁结构进行有限元分析
%单元类型：线性三角形单元
%单元数目：40X30
%单元划分思路如下：                                             
%                                                |  5000N
%                                                V
%                      /|1---------- 3 ----------5                                                                                                             
%                      /||  .        |  .        |                                                                                                              
%                      /||     .     |     .     |          
%                      /||        .  |        .  |          
%                      /|2---------- 4 ----------6 
%                                                |
%                                                v  5000N
%
%边界条件：
%左端固定端约束，右端上下节点各受到F=5000N的作用，方向如上图所示。
%%%%%%%%%
clear all;
clc;
% INITIALIZE
nelx=40;                                                                   %Number of elements in X direction
nely=30;                                                                   %Number of elements in Y direction
ly=1;                                                                      %y-direction length
lx=2;                                                                      %X-direction length
E=1*10^7;                                                                  %Young's moduli 
nu=1/3;                                                                    %Poisson's ratio
t=0.1;                                                                     %Thickness

% FE-ANALYSIS       
x0=[];
for i=1:nelx+1
    for j=1:nely+1
          x0=[x0; (i-1)*lx/nelx, -(j-1)*ly/nely];
    end
end
nodes=[];
for i=1:nelx
    for j=1:nely
        nodes=[nodes; (nely+1)*(i-1)+j, (nely+1)*(i-1)+j+1, (nely+1)*i+j+1];
        nodes=[nodes; (nely+1)*(i-1)+j, (nely+1)*i+j+1, (nely+1)*i+j];
    end
end

% Stiffness Matrix Calculation and Assembly
KK = zeros(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
for iel=1:2*nelx*nely       

    nd(1)=nodes(iel,1);         
    nd(2)=nodes(iel,2);         
    nd(3)=nodes(iel,3);        
    
    x1=x0(nd(1),1); y1=x0(nd(1),2);     
    x2=x0(nd(2),1); y2=x0(nd(2),2);     
    x3=x0(nd(3),1); y3=x0(nd(3),2); 

    k=Stiffness(E,nu,t,x1,y1,x2,y2,x3,y3);
    KK=Assembly(KK,k,nd(1),nd(2),nd(3));
end


% DEFINE LOADS AND SUPPORTS 
F(2*(nely+1)*nelx+2,1) = -5000;
F(2*(nely+1)*(nelx+1),1) = -5000;
fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

% SOLVING
U(freedofs,:) = KK(freedofs,freedofs) \F(freedofs,:);
U(fixeddofs,:)= 0;

%Stress calculation for per elements
stress=[];
for iel=1:2*nelx*nely       

    nd(1)=nodes(iel,1);         
    nd(2)=nodes(iel,2);         
    nd(3)=nodes(iel,3);        
    
    x1=x0(nd(1),1); y1=x0(nd(1),2);     
    x2=x0(nd(2),1); y2=x0(nd(2),2);     
    x3=x0(nd(3),1); y3=x0(nd(3),2); 
    
    u=[U(2*nd(1)-1); U(2*nd(1)); U(2*nd(2)-1); U(2*nd(2)); U(2*nd(3)-1); U(2*nd(3))];
    
    stress=[stress; Stress(E,nu,x1,y1,x2,y2,x3,y3,u)'];
end
    

% ELEMENT STIFFNESS MATRIX 
function k=Stiffness(E,NU,t,xi,yi,xj,yj,xm,ym)
%该函数计算单元的刚度矩阵
%输入弹性模量E，泊松比NU，厚度t
%输入三个节点i、j、m的坐标xi,yi,xj,yj,xm,ym
%输出单元刚度矩阵k(6X6)
%---------------------------------------------------------------
A = (xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2;
betai = yj-ym;
betaj = ym-yi;
betam = yi-yj;
gammai = xm-xj;
gammaj = xi-xm;
gammam = xj-xi;
B = [betai 0 betaj 0 betam 0 ;
0 gammai 0 gammaj 0 gammam ;
gammai betai gammaj betaj gammam betam]/(2*A);
D = (E/(1-NU*NU))*[1 NU 0 ; NU 1 0 ; 0 0 (1-NU)/2];
k= t*A*B'*D*B;
end

%整体刚度矩阵组装
function z = Assembly(KK,k,i,j,m)
%该函数进行单元刚度矩阵的组装
%输入单元刚度矩阵k
%输入单元的节点编号I、j、m
%输出整体刚度矩阵KK
%---------------------------------------------------------------
DOF(1)=2*i-1;
DOF(2)=2*i;
DOF(3)=2*j-1;
DOF(4)=2*j;
DOF(5)=2*m-1;
DOF(6)=2*m;
for n1=1:6
for n2=1:6
KK(DOF(n1),DOF(n2))= KK(DOF(n1),DOF(n2))+k(n1,n2);
end
end
z=KK;
end

%单元应力
function stress=Stress(E,NU,xi,yi,xj,yj,xm,ym,u)
%该函数计算单元的应力
%输入弹性模量E，泊松比NU，厚度t
%输入三个节点i、j、m的坐标xi,yi,xj,yj,xm,ym
%输出单元的应力stress(3X1)，由于它为常应力单元，则单元的应力分量为Sx,Sy,Sxy
%---------------------------------------------------------------
A = (xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2;
betai = yj-ym;
betaj = ym-yi;
betam = yi-yj;
gammai = xm-xj;
gammaj = xi-xm;
gammam = xj-xi;
B = [betai 0 betaj 0 betam 0 ;
0 gammai 0 gammaj 0 gammam ;
gammai betai gammaj betaj gammam betam]/(2*A);
D = (E/(1-NU*NU))*[1 NU 0 ; NU 1 0 ; 0 0 (1-NU)/2];
stress = D*B*u;
end

