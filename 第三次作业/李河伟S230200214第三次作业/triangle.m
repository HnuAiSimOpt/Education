%% ѧ�ţ�230200214
%���������ΰ
%���ø˵�Ԫ����Ϊ60mm����Ϊ20mm�����Ϊ5mm�����ɱ�Ϊ0.3�������غ�Ϊ
%1000N������������˵��е㴦����������,�Ҷ˹̶���
clc;
%% �������
global nnd nel nne nodof eldof sdof n geom connec dee nf Nodal_loads
format long
%% �ڵ�&��Ԫ��
nnd = 21 ;               	% �ڵ��� 				
nel = 24 ;                	% ��Ԫ��			
nne = 3 ;                	% ��Ԫ�ڵ���			
nodof =2;                	% �ڵ����ɶ�
eldof = nne*nodof;       	% ��Ԫ���ɶ�
sdof = nnd*nodof;
%% �ڵ�����
% geom(nnd,2)�洢�ڵ��x,y����
geom=[];
for i = 1 : 7
    geom=[geom; 10*(i-1), -10];
    geom=[geom; 10*(i-1), 0];
    geom=[geom; 10*(i-1), 10];
end

%% �ڵ������ϵ
%connec(nel,2)�洢�ڵ������ϵ�����ڵ������ϵ����
connec=[];
for i = 1 : 6
    for j = 1 : 2
        connec=[connec; 3*i+j-3, 3*i+j, 3*i+j-2];
        connec=[connec; 3*i+j, 3*i+j+1, 3*i+j-2];
    end
end

%% ��Ԫ���Ϻͼ��β���
E = 200000.;     % E��Mpa��
vu = 0.3;        % ���ɱ�
thick = 5.;      % ������ȣ�mm��

%% �γ�ƽ��Ӧ���ĵ��Ծ���
dee = E/(1-vu*vu).*[1, vu, 0; vu, 1, 0; 0,0,(1-vu)*0.5];

%% �߽�����
%nf(nnd,nodof)����߽�������0�������ɶȱ�Լ����1����ţ�����û�С�
%% 1.����ÿ���ڵ�ÿ�����ɶȶ�δ��Լ����
nf = ones(nnd, nodof);    % ������ nf ��ʼ��Ϊ 1
%% 2.��ʵ����Լ�������ɶ�����Ϊ0
nf(19,1) = 0; nf(19,2) = 0;  % �ڵ�19��20��21Ϊ�̶�Լ�������ɶ�0
nf(20,1) = 0; nf(20,2) = 0;  
nf(21,1) = 0; nf(21,2) = 0;  
%%3.��δ��Լ���Ľڵ���б��
n=0;
for i=1:nnd
    for j=1:nodof
        if nf(i,j) ~= 0 
            n=n+1;
           nf(i,j)=n;
        end
    end
end
 
%% �ڵ��غ�
%load���󴢴����
Nodal_loads= zeros(nnd, 2);
Nodal_loads(2,1) = 0.; Nodal_loads(2,2) = -1000.;    % Node 2

%% 
fg = zeros(n,1); %������ʸ����װ
for i = 1 : nnd
    if nf(i,1) ~= 0
        fg(nf(i,1))= Nodal_loads(i,1);
    end
    if nf(i,2) ~= 0
       fg(nf(i,2))= Nodal_loads(i,2);
    end
end

%�նȾ������װ��
kk = zeros(sdof,sdof);
nd=zeros(3,1);
for i = 1:nel
    nd(1) = connec(i,1);
    nd(2) = connec(i,2);
    nd(3) = connec(i,3);
    x1 = geom(nd(1),1);y1 = geom(nd(1),2);
    x2 = geom(nd(2),1);y2 = geom(nd(2),2);
    x3 = geom(nd(3),1);y3 = geom(nd(3),2);
    A = 0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2);
    A2 = A*2;
    dhdx=(1/A2)*[(y2-y3), (y3-y1), (y1-y2)];
    dhdy=(1/A2)*[(x3-x2), (x1-x3), (x2-x1)];
    %B����
    kinmtx2=fekine2d(nne,dhdx,dhdy);
    k=A*thick*kinmtx2'*dee*kinmtx2;
    index=feeldof(nd,nne ,nodof);
    for j1 = 1:length(index)
        ii = index(j1);
        for j2 = 1:length(index)
            jj = index(j2);
            kk(ii,jj)=kk(ii,jj)+k(j1,j2);
        end
    end
end  
%% ����
kk1=kk;
%--------------------------------------------------------------------------
%   force vector
%--------------------------------------------------------------------------
ff=fg;
%------------------------
%  apply boundary condition
%------------------------
[kk,ff]=feaplyc3(kk,ff,nf,nnd,nodof);
%[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);

%-------------------------
%  solve the matrix equation
%-------------------------

%disp=kk\ff;
[LL,UU]=lu(kk);
utemp=LL\ff;
disp=UU\utemp;
EU=0.5*disp'*kk1*disp;
%---------------------------------------
%  element stress computation
%---------------------------------------
energy=0;
stress=zeros(nel,3);		
strain=zeros(nel,3);
for ielp=1:nne           % loop for the total number of elements
     nd(1)=connec(ielp,1); % 1st connected node for (iel)-th element
     nd(2)=connec(ielp,2); % 2nd connected node for (iel)-th element
     nd(3)=connec(ielp,3); % 3rd connected node for (iel)-th element

    x1=geom(nd(1),1); y1=geom(nd(1),2);% coord values of 1st node
    x2=geom(nd(2),1); y2=geom(nd(2),2);% coord values of 2nd node
    x3=geom(nd(3),1); y3=geom(nd(3),2);% coord values of 3rd node
    
    xcentre=(x1+x2+x3)/3; ycentre=(y1+y2+y3)/3;
    
    index=feeldof(nd,nne,nodof);% extract system dofs associated with element
    for i=1:eldof
        eldisp(i)=disp(index(i));
    end
    area=0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2);  % area of triangule
    area2=area*2;
    dhdx=(1/area2)*[(y2-y3) (y3-y1) (y1-y2)];  % derivatives w.r.t. x-axis
    dhdy=(1/area2)*[(x3-x2) (x1-x3) (x2-x1)];  % derivatives w.r.t. y-axis
    kinmtx2=fekine2d(nne,dhdx,dhdy);          % compute kinematic matrix
    estrain=kinmtx2*eldisp';             % compute strains
    estress=dee*estrain;             % compute stresses
    for i=1:3
        strain(ielp,i)=estrain(i);          % store for each element
        stress(ielp,i)=estress(i);          % store for each element          
    end
    energy=energy+0.5*estrain'*dee*estrain*area;
end
nnode=nnd;
nodes=connec;
neigh_node = cell(nnode,1);
indneigh=zeros(1,nnode);
for i=1:nel
    for j=1:3
        indneigh(nodes(i,j))=indneigh(nodes(i,j))+1;
        neigh_node{nodes(i,j)}(indneigh(nodes(i,j)))=i;
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

fid_out=fopen('result_beam01.plt','w');

fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
fprintf(fid_out,'VARIABLES="x" "y" "u" "v" "sigax"  "sigmay" "sigmaxy"\n');
fprintf(fid_out,'ZONE T="flow-field", N= %8d,E=%8d,ET=TRIANGLE, F=FEPOINT\n',nnode,nel);
for i=1:nnode
       fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',geom(i,1),geom(i,2),disp(2*i-1),disp(2*i),stress_node(1,i),stress_node(2,i),stress_node(3,i));
end
for i=1:nel
      fprintf(fid_out,'%8d%8d%8d\n',nodes(i,1),nodes(i,2),nodes(i,3));
end

function [kk,ff]=feaplyc3(kk,ff,nf,nnd,nodof)
%----------------------------------------------------------
%  Purpose:
%     Apply constraints to matrix equation [kk]{x}={ff}
%  Variable Description:
%     kk - system matrix before applying constraints 
%     ff - system vector before applying constraints
%     bcdof - a vector containging constrained d.o.f
%     bcval - a vector containing contained value 
%-----------------------------------------------------------
 sdof=size(kk); 
 for i=1:nnd
     for j1=1:nodof
         if nf(i,j1)==0
             c=(i-1)*nodof+j1;
             for  j2=1:sdof
                kk(c,j2)=0;
             end
             kk(c,c)=1;
             ff(c)=nf(i,j1);
         end
     end
 end
end

function [index] = feeldof(nd,nnel,ndof)
%----------------------------------------------------------
%  Purpose:
%     Compute system dofs associated with each element 
%  Variable Description:
%     index - system dof vector associated with element "iel��
%     nnel - number of nodes per element
%     ndof - number of dofs per node 
%-----------------------------------------------------------
 k=0;
   for i=1:nnel
        start = (nd(i)-1)*ndof;
        for j=1:ndof
              k=k+1;
              index(k)=start+j;
       end
   end
end

function [kk]=feasmb_2(kk,k,index)
%----------------------------------------------------------
%  Purpose:
%     Assembly of element matrices into the system matrix
%  Variable Description:
%     kk - system matrix
%     k  - element matri
%     index - d.o.f. vector associated with an element
%-----------------------------------------------------------
edof = length(index);
for i=1:edof
    ii=index(i);
    for j=1:edof
        jj=index(j);
        kk(ii,jj)=kk(ii,jj)+k(i,j);
    end
end
end

function [kinmtx2]=fekine2d(nnel,dhdx,dhdy)
%------------------------------------------------------------------------
%  Purpose:
%     determine the kinematic equation between strains and displacements
%     for two-dimensional solids
%  Variable Description:
%     nnel - number of nodes per element
%     dhdx - derivatives of shape functions with respect to x   
%     dhdy - derivatives of shape functions with respect to y
%------------------------------------------------------------------------

for i=1:nnel
       i1=2*i-1;  
       i2=i1+1;
       kinmtx2(1,i1)=dhdx(i);
       kinmtx2(2,i2)=dhdy(i);
       kinmtx2(3,i1)=dhdy(i);
       kinmtx2(3,i2)=dhdx(i);
 end
end

