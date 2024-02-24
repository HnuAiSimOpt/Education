%-------------------------------------------------------------------------
%Problem: Cantilever beam : quadrilateral element 
% Using FEM-Q4 element
% Developed by Cui Xiangyang, Hunan University 2008/08
%--------------------------------------------------------------------------
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
clear
format short
first_time=cputime; 
%---------------------------------
%input data for control parameters
lengthx=4;              %length of x-axis side of problem
lengthy=4;              %length of y-axis side of problem

lx=48;                  % number of element in x-axis
ly=24;                   % number of element in y-axis

nel=lx*ly;               % number of element

nnel=4;                 %number of nodes per element
ndof=2;                 %number of dofs per node
nnode=(lx+1)*(ly+1);    %total number of nodes in system    
sdof=nnode*ndof;        %total system dofs
edof=nnel*ndof;         %degrees of freedom per element

emodule=1.0;            %elastic modulus
poisson=0.0;            %Poisson's ratio

nglx=2; ngly=2;     %2x2 Gauss-Legendre quadrature
nglxy=nglx*ngly;    %number of sampling points per element

fload=-1;         % the total load 

%--------------------------------------------
%input data for nodal coordinate values
%gcoord(i,j) where i->node no.  and j->x or y
%--------------------------------------------

x0=[];
for i=1:lx+1
    for j=1:ly+1
        x0=[x0; (i-1)*lengthx/lx  (j-1)*lengthy/ly];%-0.5*lengthy*(1+(lx+1-i)/lx)*(1-(j-1)/ly)
    end
end

%-----------------------------------------------------
%input data for nodal connectivity for each element
%nodes(i,j) where i->element no. and j->connected nodes
%------------------------------------------------------

nodes=[];
for i=1:lx
    for j=1:ly
        nodes=[nodes; (ly+1)*(i-1)+j (ly+1)*i+j (ly+1)*i+j+1 (ly+1)*(i-1)+j+1;];
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
   for ie=1:nel
        for j=1:nnel+1
            j1=mod(j-1,nnel)+1;
            xp(j)=x0(nodes(ie,j1),1);
            yp(j)=x0(nodes(ie,j1),2);
        end
        plot(xp,yp,'-')
   end
%   pause           
end

%----------------------------------
%input data for boundary conditions
%----------------------------------

bcdof=[];
bcval=[];
for i=1:ly+1
        bcdof=[bcdof 1+2*(i-1) 2+2*(i-1)];
        bcval=[bcval 0 0];
end
        
%--------------------------------------
%initialization of matrices and vectors
%--------------------------------------
fid=fopen('testquadstress.txt','w');
for i = 1 : nnode
    if abs(x0(i, 1)) < 1e-3)
        fprintf(fid,'node %d  coords 2 %e %e boundary\r\n',i,x0(i,1), x0(i,2));
    else 
        fprintf(fid,'node %d  coords 2 %e %e \r\n',i,x0(i,1), x0(i,2)); 
    end
end
for i = 1 : nel
    fprintf(fid,'quadplanestress %d nodes 4 %d %d %d %d\r\n', i, nodes(i,1), nodes(i,2), nodes(i,3), nodes(i,4));
end
fclose(fid);

ff=sparse(sdof,1);			%system force vector
k=sparse(edof,edof);		%initialization of element matrix
kk=zeros(sdof,sdof);		%system matrix
disp=sparse(sdof,1);		%system displacement vector
eldisp=sparse(edof,1);		%element displacement vector
stress=zeros(nel,4,3);		%matrix containing stress components
strain=sparse(nel,4,3);		%matrix containing strain components
index=sparse(edof,1);		%index vector
kinmtx=sparse(3,edof);		%kinematic matrix
matmtx=sparse(3,3);			%constitutive matrix

%-------------------------------------------------
%compute element matrices and vectors and assemble
%-------------------------------------------------

[point2,weight2]=feglqd2(nglx,ngly);    %sampling points & weights 高斯点的坐标和权重
matmtx=fematiso(1,emodule,poisson);     %constitutive matrice
for iel=1:nel       %loop for the total number of element 
    for i=1:nnel    %取出单元中每个节点的坐标值
        nd(i)=nodes(iel,i);             %extract nodes for (iel)-th element
        xcoord(i)=x0(nd(i),1);      %extract x value of the node
        ycoord(i)=x0(nd(i),2);      %extract y value of the node
    end
    s=(ycoord(3)-ycoord(2)+ycoord(4)-ycoord(1))*xcoord(2)/2;
    [angleA,angleB,angleC,angleD]=cal_angle(xcoord,ycoord); %4个角的角度
    if angleA >=pi | angleA >=pi | angleA >=pi | angleA >=pi 
        temp='quadrilateral element is concave - FEM can not solve'
        break;
    end
    k=sparse(edof,edof);			%initialization of element matrix
     
    %---------------------
    %numerical integration
    %---------------------
    for intx=1:nglx
        x=point2(intx,1);           %sampling point in x-axis
        wtx=weight2(intx,1);        %weight in x-axis
        for inty=1:ngly
            y=point2(inty,2);       %sampling point in y-axis
            wty=weight2(inty,2);    %weight in y-axis
            
            [shape,dhdr,dhds]=feisoq4(x,y);     %compute shape functions and derivatives at sampling point
            
            jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  %compute Jacobian
            
            detjacob=det(jacob2);                           %determinant of Jacobian           
           %% mm=mm+detjacob
            invjacob=inv(jacob2);                           %inverse of Jacobian matrix
            
            [dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob);  %derivatives w.r.t physical coordinate
            
            kinmtx2=fekine2d(nnel,dhdx,dhdy);               %kinematic matrice
   
            k=k+kinmtx2'*matmtx*kinmtx2*wtx*wty*detjacob;   %element stiffness matrice
            
        end
    end
    
    index=feeldof(nd,nnel,ndof);        %extract system dofs for the element
    
   	kk=feasmb_2(kk,k,index);            %assemble element matrics
end

kk1=kk;
%--------------------------------------------------------------------------
%force vector
%--------------------------------------------------------------------------
ff(sdof,1)=fload;
%------------------------
%apply boundary condition
%------------------------
[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);

%-------------------------
%solve the matrix equation
%-------------------------

%disp=kk\ff;
[LL UU]=lu(kk);
utemp=LL\ff;
disp=UU\utemp;

solve_time = cputime-first_time;

EU=0.5*disp'*kk1*disp;
U_fem=disp;
%-------------------------------------------------
%  element stress computation
%-------------------------------------------------
energy=0;
mk = zeros(4, 4);
count = 1;
for intx=1:nglx
        x=point2(intx,1);           %sampling point in x-axis
        for inty=1:ngly
            y=point2(inty,2);       %sampling point in y-axis
            [shape,dhdr,dhds]=feisoq4(x,y);
            for i = 1 : 4 
                mk(count, i) = shape(i);
            end
            count = count + 1;
        end
end
stab2 = inv(mk);
stab=[1.866 -0.5 0.134 -0.5;
      -0.5 1.866 -0.5 0.134;
      0.134 -0.5 1.866 -0.5;
      -0.5 0.134 -0.5 1.866];
  
for iel=1:nel       %loop for the total number of element
    for i=1:nnel
        nd(i)=nodes(iel,i);             %extract nodes for (iel)-th element
        xcoord(i)=x0(nd(i),1);      %extract x value of the node
        ycoord(i)=x0(nd(i),2);      %extract y value of the node
    end
    [angleA,angleB,angleC,angleD]=cal_angle(xcoord,ycoord);
    if angleA >=pi | angleA >=pi | angleA >=pi | angleA >=pi
        temp='quadrilateral element is concave - FEM can not solve'
        break; 
    end    
    k=sparse(edof,edof);			%initialization of element matrix   
    
    %---------------------
    %numerical integration
    %---------------------
    pk=0;
    gstress=zeros(3,4);
    for intx=1:nglx
        x=point2(intx,1);           %sampling point in x-axis
        wtx=weight2(intx,1);        %weight in x-axis
        for inty=1:ngly
            y=point2(inty,2);       %sampling point in y-axis
            wty=weight2(inty,2);    %weight in y-axis
            
            [shape,dhdr,dhds]=feisoq4(x,y);     %compute shape functions and derivatives at sampling point
            
            jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  %compute Jacobian
            
            detjacob=det(jacob2);                           %determinant of Jacobian           
            
            invjacob=inv(jacob2);                           %inverse of Jacobian matrix
            
            [dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob);  %derivatives w.r.t physical coordinate
            
            kinmtx2=fekine2d(nnel,dhdx,dhdy);               %kinematic matrice
            
            index=feeldof(nd,nnel,ndof);  %extract system dofs for the element
            
            for i=1:edof
                eldisp(i,1)=disp(index(i));
            end            
            estrain=kinmtx2*eldisp;             % compute strains
            estress=matmtx*estrain;             % compute stresses           
            energy=energy+0.5*estrain'*matmtx*estrain*wtx*wty*detjacob;
            pk=pk+1;
            gstress(:,pk)=estress;
        end
    end
    for i=1:3
        for j=1:4
            for n=1:4
                stress(iel,j,i)= stress(iel,j,i)+stab2(j,n)*gstress(i,n);
            end
        end
    end
    
end

neigh_node = cell(nnode,1);
neigh_node_ind = cell(nnode,1);
indneigh=zeros(1,nnode);
for i=1:nel
    for j=1:4
        indneigh(nodes(i,j))=indneigh(nodes(i,j))+1;
        neigh_node{nodes(i,j)}(indneigh(nodes(i,j)))=i;
        neigh_node_ind{nodes(i,j)}(indneigh(nodes(i,j)))=j;
    end
end

stress_node=zeros(3,nnode);	
for inode=1:nnode
    numel= indneigh(inode);
    for i=1:numel
        ind_nel= neigh_node{inode}(i);
        ind_nod=neigh_node_ind{inode}(i);
        for j=1:3
            stress_node(j,inode)=stress_node(j,inode)+stress(ind_nel,ind_nod,j);
        end
    end
    stress_node(:,inode)=stress_node(:,inode)/numel;
end
%-------------------------------------------------------------
% Data output
%-------------------------------------------------------------
fid_out=fopen('result_beamQ4_2.plt','w');
fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
fprintf(fid_out,'VARIABLES="x" "y" "u" "v" "sigax"  "sigmay" "sigmaxy"\n');
fprintf(fid_out,'ZONE T="flow-field", N= %8d,E=%8d,ET=QUADRILATERAL, F=FEPOINT\n',nnode,nel);
for i=1:nnode
    fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',x0(i,1),x0(i,2),disp(2*i-1),disp(2*i),stress_node(1,i),stress_node(2,i),stress_node(3,i));
end
for i=1:nel
    fprintf(fid_out,'%8d%8d%8d%8d\n',nodes(i,1),nodes(i,2),nodes(i,3),nodes(i,4));
end