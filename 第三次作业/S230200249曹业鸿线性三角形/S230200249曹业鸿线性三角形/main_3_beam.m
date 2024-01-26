%------------------------------------
%Problem: Cantilever beam , triangle element
% undeted by Cui Xiangyang
%Variable descriptions
% k = element stiffness matrix
% f = element force vector
% kk = system stiffness matrix
% ff = system force vector
% disp = system nodal displacement vector
% eldisp = element nodal displacement vector
% stress = matrix containing stresses
% strain = matrix containing strains
% x0 = coordinate values of each node
% nodes = nodal connectivity of each element
% index = a vector containing systems dofs associated with element
% bcdof = a vector containing dofs associated with boundary conditions
% bcval = a vector containing boundary condition values associated with the dofs in bcdof
% ---------------------------------------------------------------------------------------

clear all

first_time=cputime; 

format long

%---------------------------------
%input data for control parameters
%---------------------------------

lengthx=4;              %length of x-axis side of problem
lengthy=2;              %length of y-axis side of problem

lx=48;                  % number of element in x-axis
ly=24;                   % number of element in y-axis
nel=2*lx*ly;            % number of element

nnel=3;                 %number of nodes per element
ndof=2;                 %number of dofs per node
nnode=(lx+1)*(ly+1);    %total number of nodes in system    
sdof=nnode*ndof;        %total system dofs
edof=nnel*ndof;         %degrees of freedom per element
emodule=1.0;            %elastic modulus
poisson=0.0;            %Poisson's ratio

fload=-1;             % the total load 

%--------------------------------------------
%input data for nodal coordinate values
%gcoord(i,j) where i->node no.  and j->x or y
%--------------------------------------------

x0=[];
for i=1:lx+1
    for j=1:ly+1
        x0=[x0; (i-1)*lengthx/lx  -0.5*lengthy*(1+(lx+1-i)/lx)*(1-(j-1)/ly)];
    end
end

%-----------------------------------------------------
%input data for nodal connectivity for each element
%nodes(i,j) where i->element no. and j->connected nodes
%------------------------------------------------------

nodes=[];
for i=1:lx
    for j=1:ly
        nodes=[nodes; (ly+1)*(i-1)+j (ly+1)*i+j (ly+1)*(i-1)+j+1;];
        nodes=[nodes; (ly+1)*i+j (ly+1)*i+j+1 (ly+1)*(i-1)+j+1;];
    end
end

%**************************************************************************
% Output the mesh
iplot=1;
if iplot==1
   figure(99)
   hold on
   axis off
   axis equal
   for ie=1:nel  %����Ԫ����ӡ
        for j=1:nnel+1 %j=1 2 3 4
            j1=mod(j-1,nnel)+1; %j1=1 2 3 1  �պý�һ��������������
            xp(j)=x0(nodes(ie,j1),1);
            yp(j)=x0(nodes(ie,j1),2);
        end
        plot(xp,yp,'-')
   end
%   pause           
end
%--------------------------------------------------------------------------
%----------------------------------
%input data for boundary conditions
%----------------------------------

bcdof=[]; %�洢���Ǳ�Լ���Ľڵ�ı��
bcval=[]; %�洢����Լ����ֵ�Ĵ�С
for i=1:ly+1
        bcdof=[bcdof 1+2*(i-1) 2+2*(i-1)];
        bcval=[bcval 0 0];
end
        
%--------------------------------------
%initialization of matrices and vectors  ����һϵ�о���
%--------------------------------------

ff=sparse(sdof,1);			%system force vector
k=sparse(edof,edof);			%initialization of element matrix
kk=sparse(sdof,sdof);		%system matrix
disp=sparse(sdof,1);			%system displacement vector
eldisp=sparse(edof,1);		%element displacement vector
stress=zeros(nel,3);		%matrix containing stress components
strain=zeros(nel,3);		%matrix containing strain components
index=sparse(edof,1);		%index vector
kinmtx=sparse(3,edof);		%kinematic matrix
matmtx=sparse(3,3);			%constitutive matrix

%-------------------------------------------------
%compute element matrices and vectors and assemble
%-------------------------------------------------

matmtx=fematiso(1,emodule,poisson);     %constitutive matrice������

for iel=1:nel       %loop for the total number of element ��Ԫ���

    nd(1)=nodes(iel,1);         %1st connected node for (iel)-th element
    nd(2)=nodes(iel,2);         %2nd connected node for (iel)-th element
    nd(3)=nodes(iel,3);         %3rd connected node for (iel)-th element ��Ԫ�Ľڵ���
    
    x1=x0(nd(1),1); y1=x0(nd(1),2);     %coord values of 1st node
    x2=x0(nd(2),1); y2=x0(nd(2),2);     %coord values of 2nd node
    x3=x0(nd(3),1); y3=x0(nd(3),2);     %coord values of 3rd node   �ڵ��Ӧ��������Ϣ
    
    index=feeldof(nd,nnel,ndof);        %extract system dofs for the element ��Ԫ�ն��е��������ܸնȾ����е���������
    
    %---------------------------------------
    %find the derivatives of shape functions
    %---------------------------------------
    
   	area=0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2); %area of triangula
	area2=area*2;
	
    dhdx=(1/area2)*[(y2-y3) (y3-y1) (y1-y2)];       %derivatives w.r.t x //bi bj bm
    dhdy=(1/area2)*[(x3-x2) (x1-x3) (x2-x1)];       %derivatives w.r.t y //ci cj cm
   
 	kinmtx2=fekine2d(nnel,dhdx,dhdy);               %kinematic matrice
   
   	k=kinmtx2'*matmtx*kinmtx2*area;                 %element stiffness matrice
     
   	kk=feasmb_2(kk,k,index);                        %assemble element matrics
end

kk1=kk;
%--------------------------------------------------------------------------
%force vector
%--------------------------------------------------------------------------
ff(sdof,1)=fload;  %�������ؾ���ÿ���ڵ�ĸ���Ϊ-1
%------------------------
%apply boundary condition
%------------------------
[kk,ff]=feaplyc2(kk,ff,bcdof,bcval); %���ܸնȾ���͸��ؾ������Լ���� ����ߵ�һ�е�Ԫλ��Ϊ0

%-------------------------
%solve the matrix equation
%-------------------------

%disp=kk\ff;
[LL UU]=lu(kk);       %��kk�ֽ�Ϊһ�������Ǻ������Ǿ���ĳ˻���������и�˹��Ԫ��
utemp=LL\ff;
disp=UU\utemp;         %��������ڵ��λ��

solve_time = cputime-first_time;

EU=0.5*disp'*kk1*disp     %����
U_fem=disp;
%save d:\Make_program\FEM\FEM_2D_standard\Cantilever_beam\U_FEM U_fem;
%---------------------------------------
%  element stress computation
%---------------------------------------
energy=0;
energy0=0;
denergy=0;
for ielp=1:nel           % loop for the total number of elements
    nd(1)=nodes(ielp,1); % 1st connected node for (iel)-th element
    nd(2)=nodes(ielp,2); % 2nd connected node for (iel)-th element
    nd(3)=nodes(ielp,3); % 3rd connected node for (iel)-th element
    %��ie����Ԫ�Ľڵ���

    x1=x0(nd(1),1); y1=x0(nd(1),2);% coord values of 1st node
    x2=x0(nd(2),1); y2=x0(nd(2),2);% coord values of 2nd node
    x3=x0(nd(3),1); y3=x0(nd(3),2);% coord values of 3rd node
    %�����Ԫ�ڵ��������Ϣ
    xcentre=(x1+x2+x3)/3; ycentre=(y1+y2+y3)/3;
    %��Ԫ����(���ģ�
    index=feeldof(nd,nnel,ndof);% extract system dofs associated with element
    %��Ԫ�նȾ����е�Ԫ�����ܸնȾ����е�λ��
    %-------------------------------------------------------
    %  extract element displacement vector  %��ȡ��Ԫλ������
    %-------------------------------------------------------

    for i=1:edof
        eldisp(i)=disp(index(i));
    end          %����ie���ڵ��λ��ֵ(u v)����eldisp

    area=0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2);  % area of triangule
    area2=area*2;
    dhdx=(1/area2)*[(y2-y3) (y3-y1) (y1-y2)];  % derivatives w.r.t. x-axis
    dhdy=(1/area2)*[(x3-x2) (x1-x3) (x2-x1)];  % derivatives w.r.t. y-axis

    kinmtx2=fekine2d(nnel,dhdx,dhdy);          % compute kinematic matrix

    estrain=kinmtx2*eldisp;             % compute strains
    estress=matmtx*estrain;             % compute stresses

    for i=1:3
        strain(ielp,i)=estrain(i);          % store for each element
        stress(ielp,i)=estress(i);          % store for each element          
    end
    
    energy=energy+0.5*estrain'*matmtx*estrain*area;
end

neigh_node = cell(nnode,1);%����nnode*1��Ԫ������
indneigh=zeros(1,nnode);
for i=1:nel
    for j=1:3
        indneigh(nodes(i,j))=indneigh(nodes(i,j))+1; %����ڵ�ĵ�Ԫ��
        neigh_node{nodes(i,j)}(indneigh(nodes(i,j)))=i;%����ڵ�ĵ�Ԫ���
    end
end
stress_node=zeros(3,nnode);	
for inode=1:nnode
    numel= indneigh(inode);
    for i=1:numel
        ind_nel= neigh_node{inode}(i);
        for j=1:3
            stress_node(j,inode)=stress_node(j,inode)+stress(ind_nel,j);
        end
    end
    stress_node(:,inode)=stress_node(:,inode)/numel;
end
%-------------------------------------------------------------
% Data output
%-------------------------------------------------------------
fid_out=fopen('result_beam01.plt','w');
fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
fprintf(fid_out,'VARIABLES="x" "y" "u" "v" "sigax"  "sigmay" "sigmaxy"\n');
fprintf(fid_out,'ZONE T="flow-field", N= %8d,E=%8d,ET=TRIANGLE, F=FEPOINT\n',nnode,nel);
for i=1:nnode
    fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',...
        x0(i,1),x0(i,2),disp(2*i-1),disp(2*i),stress_node(1,i),stress_node(2,i),stress_node(3,i));
end
for i=1:nel
    fprintf(fid_out,'%8d%8d%8d\n',nodes(i,1),nodes(i,2),nodes(i,3));
end
