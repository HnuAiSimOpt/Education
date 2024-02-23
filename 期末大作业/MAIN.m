% 二阶四面体单元
% 模型文件C3D10Test.txt
%Variable descriptions
% k = element stiffness matrix
% kk = system stiffness matrix
% ff = system force vector
% disp = system nodal displacement vector
% eldisp = element nodal displacement vector
% xy = coordinate values of each node
% nodes = nodal connectivity of each element
% index = a vector containing systems dofs associated with element
% bcdof = a vector containing dofs associated with boundary conditions
% bcval = a vector containing boundary condition values associated with the dofs in bcdof 

clear all;
clc;
format long

%-----------------------------------------------------
%input data for nodal connectivity for each element
%nodes(i,j) where i->element no. and j->connected nodes
%------------------------------------------------------
mesh=dlmread('C3D10Test.txt',',');
nnel=10;                 %number of nodes per element
ndof=3;                 %number of dofs per node
for i=1:length(mesh(:,1))
    if mesh(i,1)==1&&mesh(i,5)~=0
        flag=i;
        break;
    end
end
xy=mesh(1:flag-1,1:ndof+1);
nodes=mesh(flag:length(mesh(:,1)),1:nnel+1);

% PAPERING WORK
nel = length(nodes(:,1));               % number of elements of system.
nnode = length(xy(:,1));             % number of nodes of system.
sdof = nnode * ndof;                       % number of DOFs of system.
edof = ndof * nnel;                      % number of DOFs per element.
ff = sparse(sdof,1);                    % force vector of system.
kk = sparse(sdof,sdof);                  % system stiffness matrix.  EFFECTIVE--SPARSE TYPE
disp = sparse(sdof,1);                    % displacement vector of system.
index = sparse(edof,1);                   % index vector.
%==========================================================================
elastic=2.1e6;
poisson=0.3;

%constitutive matrice
matmtx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson   0   0    0; 
   poisson  (1-poisson)   poisson   0   0    0;
   poisson  poisson  (1-poisson)    0   0    0;
   0    0    0    (1-2*poisson)/2   0    0;
   0    0    0    0    (1-2*poisson)/2   0;
   0    0    0    0    0   (1-2*poisson)/2];



%sampling points
alpha=0.58541020;
belta=0.13819660;

%-------------------------------------------------
%compute element matrices and vectors and assemble
%-------------------------------------------------
x=[];y=[];z=[];
for count1=1:nel
    k=zeros(nnel*ndof,nnel*ndof);
    cnode=nodes(count1,2:nnel+1);
    for count2=1:nnel
        x(count2)=xy(cnode(count2),2);
        y(count2)=xy(cnode(count2),3);
        z(count2)=xy(cnode(count2),4);
    end
    
    V6=[1,1,1,1;
        x(1),x(2),x(3),x(4);
        y(1),y(2),y(3),y(4);
        z(1),z(2),z(3),z(4)];
    detv6=det(V6);
    V=detv6/6;
    
    %---------------------------------------
    %find the derivatives of shape functions
    %---------------------------------------
    
    for count3=1:4 %积分点数目为4个
        
        jfd=[belta;belta;belta;belta;];
        jfd(count3)=alpha;
        L1=jfd(1);
        L2=jfd(2);
        L3=jfd(3);
        L4=jfd(4);
        
        dndl=zeros(3,10);% 形函数导数初始化
        
        dndl(1,1)=4*L1-1;
        dndl(2,2)=4*L2-1;
        dndl(3,3)=4*L3-1;
        dndl(1,4)=1-4*L4;
        dndl(2,4)=1-4*L4;
        dndl(3,4)=1-4*L4;

        dndl(1,5)=4*L2;
        dndl(2,5)=4*L1;
        
        dndl(1,7)=4*L3;
        dndl(3,7)=4*L1;
        
        dndl(1,8)=4*L4-4*L1;
        dndl(2,8)=-4*L1;
        dndl(3,8)=-4*L1;
        
        dndl(2,6)=4*L3;
        dndl(3,6)=4*L2;
        
        dndl(1,10)=-4*L3;
        dndl(2,10)=-4*L3;
        dndl(3,10)=4*L4-4*L3;
        
        dndl(1,9)=-4*L2;
        dndl(2,9)=4*L4-4*L2;
        dndl(3,9)=-4*L2;

        Jacobi=dndl*[x',y',z'];
        detJ=det(Jacobi);
        dndx=inv(Jacobi)*dndl;
        B=[];
        for count=1:nnel
            Bt=zeros(6,3);
            Bt(1,1)=dndx(1,count);
            Bt(2,2)=dndx(2,count);
            Bt(3,3)=dndx(3,count);
            Bt(4,1)=Bt(2,2);
            Bt(4,2)=Bt(1,1);
            Bt(5,2)=Bt(3,3);
            Bt(5,3)=Bt(2,2);
            Bt(6,1)=Bt(3,3);
            Bt(6,3)=Bt(1,1);
            B=[B,Bt];
        end
        k=k+B'*matmtx*B*V/4;
    end
    index=feeldof(cnode,nnel,ndof);
    kk=feasmbl(kk,k,index);
end

%--------------------------------------------------------------------------
%force vector
%--------------------------------------------------------------------------
ff(1)=100;
ff(2)=-200;

%----------------------------------
%input data for boundary conditions
%----------------------------------
BCnode=[3,4,6,8,10,11,42,43,100,101,158,159,160,161,1053,1063,1397,1398,1399,1906,1907,1909,1985,2032,2034,2035,2039,2050,2051,2052,2058,2060,2062,2070,2110,2124,2149,2246,2695];
bcdof=[];
bcval=zeros(sdof,1);
for i=1:length(BCnode)
    bcdof=[bcdof,3*BCnode(i)-2,3*BCnode(i)-1,3*BCnode(i)];
end

%------------------------
%apply boundary condition
%------------------------
[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);

%-------------------------
%solve the matrix equation
%-------------------------
fprintf('LU decomposition......\n');
[LL,UU] = lu(kk);                                       % LU分解==>LU = kk   LL*UU*disp = ff.
fprintf('Solving......\n');
Utemp = LL\ff;                                          % Utemp = UU*disp     
fprintf('Solving......\n');
disp = UU\Utemp;

%-------------------------------------------------------------
% Data output
%-------------------------------------------------------------

x=[];
y=[];
z=[];
for i=1:nnode
    x(i)=disp(3*i-2);
    y(i)=disp(3*i-1);
    z(i)=disp(3*i);
    M(i)=sqrt(x(i).^2+y(i).^2+z(i).^2);
end
fid = fopen('Results.plt','w');
fprintf(fid,'TITLE="erjiesimianti"\n');
fprintf(fid,'VARIABLES="X""Y""Z""U""V""W""M"\n');
maxnode=max([max(nodes(:,2)),max(nodes(:,3)),max(nodes(:,4)),max(nodes(:,5))]);
fprintf(fid,'ZONE T="flow-field",N=%8d,E=%8d,ET=TETRAHEDRON,F=FEPOINT\n',maxnode,nel);
for count = 1:maxnode
    fprintf(fid,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',xy(count,2),xy(count,3),xy(count,4),x(count),y(count),z(count),M(count));
end

for i=1:nel
    fprintf(fid,'%8d %8d %8d %8d\n',nodes(i,2),nodes(i,3),nodes(i,4),nodes(i,5));
end

fclose(fid);
