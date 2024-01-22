%--------------------------------------------
%2D线性三角形单元编程：
%--------------------------------------------
clear all
first_time=cputime; 
format long
%--------------------------------------------
%input data for control parameters
%--------------------------------------------
lengthx=8;              %length of x-axis side of problem
lengthy=3;              %length of y-axis side of problem

emodule=1.0;            %elastic modulus
poisson=0.0;            %Poisson's ratio

fload=-1;             % the total load 
lx=64;                  % number of element in x-axis
ly=24;                   % number of element in y-axis
nel=2*lx*ly;            % number of element
nnel=3;                 %number of nodes per element
ndof=2;                 %number of dofs per node
nnode=(lx+1)*(ly+1);    %total number of nodes in system    
sdof=nnode*ndof;        %total system dofs
edof=nnel*ndof;         %degrees of freedom per element
%-------------------------------------------------------
%input data for nodal coordinate values
%-------------------------------------------------------
x0=[];
for i=1:lx+1
    for j=1:ly+1
          x0=[x0; (i-1)*lengthx/lx      -0.5*lengthy*(1+(lx+1-i)/lx)*(1-(j-1)/ly)];
    end
end
%----------------------------------------------------------------------
%input data for nodal connectivity for each element
%nodes(i,j) where i->element no. and j->connected nodes
%----------------------------------------------------------------------
nodes=[];
for i=1:lx
    for j=1:ly
        nodes=[nodes; (ly+1)*(i-1)+j (ly+1)*i+j (ly+1)*(i-1)+j+1;];
        nodes=[nodes; (ly+1)*i+j (ly+1)*i+j+1 (ly+1)*(i-1)+j+1;];
    end
end
%----------------------------------
%input data for boundary conditions
%----------------------------------
bcdof=[];
bcval=[];
for i=1:ly+1
        bcdof=[bcdof 1+2*(i-1) 2+2*(i-1)];
        bcval=[bcval  0   0];
end
%-------------------------------------------------
%initialization of matrices and vectors
%-------------------------------------------------
ff=sparse(sdof,1);			%system force vector
k=sparse(edof,edof);		%initialization of element matrix
kk=sparse(sdof,sdof);		%system matrix
disp=sparse(sdof,1);		%system displacement vector
eldisp=sparse(edof,1);		%element displacement vector
stress=zeros(nel,3);		%matrix containing stress components
strain=zeros(nel,3);		%matrix containing strain components
index=sparse(edof,1);		%index vector
kinmtx=sparse(3,edof);		%kinematic matrix
matmtx=sparse(3,3);		%constitutive matrix
%-------------------------------------------------
%compute material matrices
%-------------------------------------------------
matmtx=fematiso(1,emodule,poisson);     %constitutive matrice
%-------------------------------------------------
%compute element matrices and vectors and assemble
%-------------------------------------------------
for iel=1:nel       %loop for the total number of element

    nd(1)=nodes(iel,1);         %1st connected node for (iel)-th element
    nd(2)=nodes(iel,2);         %2nd connected node for (iel)-th element
    nd(3)=nodes(iel,3);         %3rd connected node for (iel)-th element
    
    x1=x0(nd(1),1); y1=x0(nd(1),2);     %coord values of 1st node
    x2=x0(nd(2),1); y2=x0(nd(2),2);     %coord values of 2nd node
    x3=x0(nd(3),1); y3=x0(nd(3),2);     %coord values of 3rd node
    
    index=feeldof(nd,nnel,ndof);        %extract system dofs for the element    
    %---------------------------------------
    %find the derivatives of shape functions
    %---------------------------------------
    
   	area=0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2); %area of triangula
	area2=area*2;
	
     dhdx=(1/area2)*[(y2-y3) (y3-y1) (y1-y2)];       %derivatives w.r.t x
     dhdy=(1/area2)*[(x3-x2) (x1-x3) (x2-x1)];       %derivatives w.r.t y
   
 	kinmtx2=fekine2d(nnel,dhdx,dhdy);               %kinematic matrice   
   	k=kinmtx2'*matmtx*kinmtx2*area;              %element stiffness matrice

     
   	kk=feasmb_2(kk,k,index);                               %assemble element matrics
end
kk1=kk;
%--------------------------------------------------------------------------
%   force vector
%--------------------------------------------------------------------------
ff(sdof,1)=fload;
%------------------------
%  apply boundary condition
%------------------------
[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);

%-------------------------
%  solve the matrix equation
%-------------------------

%disp=kk\ff;
[LL UU]=lu(kk);
utemp=LL\ff;
disp=UU\utemp;
EU=0.5*disp'*kk1*disp;
%---------------------------------------
%  element stress computation
%---------------------------------------
energy=0;
for ielp=1:nel           % loop for the total number of elements
     nd(1)=nodes(ielp,1); % 1st connected node for (iel)-th element
     nd(2)=nodes(ielp,2); % 2nd connected node for (iel)-th element
     nd(3)=nodes(ielp,3); % 3rd connected node for (iel)-th element

    x1=x0(nd(1),1); y1=x0(nd(1),2);% coord values of 1st node
    x2=x0(nd(2),1); y2=x0(nd(2),2);% coord values of 2nd node
    x3=x0(nd(3),1); y3=x0(nd(3),2);% coord values of 3rd node
    
    xcentre=(x1+x2+x3)/3; ycentre=(y1+y2+y3)/3;
    
    index=feeldof(nd,nnel,ndof);% extract system dofs associated with element
    %-------------------------------------------------------
    %  extract element displacement vector
    %-------------------------------------------------------
   for i=1:edof
        eldisp(i)=disp(index(i));
    end
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

%% 输出结果
%-------------------------------------------------------------
% Data output
%-------------------------------------------------------------
fid_out=fopen('result_beam01.plt','w');

fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
fprintf(fid_out,'VARIABLES="x" "y" "u" "v" "sigax"  "sigmay" "sigmaxy"\n');
fprintf(fid_out,'ZONE T="flow-field", N= %8d,E=%8d,ET=TRIANGLE, F=FEPOINT\n',nnode,nel);
for i=1:nnode
       fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',x0(i,1),x0(i,2), ...
                    disp(2*i-1),disp(2*i),stress_node(1,i),stress_node(2,i),stress_node(3,i));
end
for i=1:nel
      fprintf(fid_out,'%8d%8d%8d\n',nodes(i,1),nodes(i,2),nodes(i,3));
end

%% definition of function
function [matmtrx]=fematiso(iopt,elastic,poisson)
%------------------------------------------------------------------------
%  Purpose:
%     determine the constitutive equation for isotropic material
%  Variable Description:
%     elastic - elastic modulus
%     poisson - Poisson's ratio   
%     iopt=1 - plane stress analysis
%     iopt=2 - plane strain analysis
%     iopt=3 - axisymmetric analysis
%     iopt=4 - three dimensional analysis
%------------------------------------------------------------------------
if  iopt==1        % plane stress
 matmtrx= elastic/(1-poisson*poisson)* ...
    [1  poisson 0; ...
     poisson  1  0; ...
     0  0  (1-poisson)/2];
elseif   iopt==2        % plane strain
     matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
       [(1-poisson)  poisson 0; 
       poisson  (1-poisson)  0;
       0  0  (1-2*poisson)/2];
elseif  iopt==3     % axisymmetry
 matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson  0; 
   poisson  (1-poisson)   poisson  0;
   poisson  poisson  (1-poisson)   0;
   0    0    0   (1-2*poisson)/2];
elseif  iopt==4        % three-dimension
 matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson   0   0    0; 
   poisson  (1-poisson)   poisson   0   0    0;
   poisson  poisson  (1-poisson)    0   0    0;
   0    0    0    (1-2*poisson)/2   0    0;
   0    0    0    0    (1-2*poisson)/2   0;
   0    0    0    0    0   (1-2*poisson)/2];
end
end

function [index]=feeldof(nd,nnel,ndof)
%----------------------------------------------------------
%  Purpose:
%     Compute system dofs associated with each element 
%  Variable Description:
%     index - system dof vector associated with element "iel“
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
       i1=(i-1)*2+1;  
       i2=i1+1;
       kinmtx2(1,i1)=dhdx(i);
       kinmtx2(2,i2)=dhdy(i);
       kinmtx2(3,i1)=dhdy(i);
       kinmtx2(3,i2)=dhdx(i);
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

function [kk,ff]=feaplyc2(kk,ff,bcdof,bcval)
%----------------------------------------------------------
%  Purpose:
%     Apply constraints to matrix equation [kk]{x}={ff}
%  Variable Description:
%     kk - system matrix before applying constraints 
%     ff - system vector before applying constraints
%     bcdof - a vector containging constrained d.o.f
%     bcval - a vector containing contained value 
%-----------------------------------------------------------
 n=length(bcdof);
 sdof=size(kk); 
 for i=1:n
      c=bcdof(i);
      for  j=1:sdof
             kk(c,j)=0;
      end
    kk(c,c)=1;
    ff(c)=bcval(i);
 end
end

