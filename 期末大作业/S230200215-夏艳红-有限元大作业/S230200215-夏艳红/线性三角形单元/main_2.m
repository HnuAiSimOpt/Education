first_time=cputime; 
format long
%--------------------------------------------
%input data for control parameters
%--------------------------------------------
lengthx=2;              %length of x-axis side of problem
lengthy=2;              %length of y-axis side of proble
h=0.2;                  %thickness of problem
emodule=2.1*10^11;            %elastic modulus
poisson=0.3;            %Poisson's ratio
fload=-12500;             % the total load
lx=16;                  % number of element in x-axis
ly=16;                   % number of element in y-axis
nel=2*0.75*lx*ly;            % number of element
nnel=3;                 %number of nodes per element
ndof=2;                 %number of dofs per node
nnode=(lx/2+1)*(ly+1)+(lx/2)*(ly/2+1);    %total number of nodes in system    
sdof=nnode*ndof;        %total system dofs
edof=nnel*ndof;         %degrees of freedom per element
x0=[];
for i=1:lx+1
    if i<=(lx/2+1)
         for j=1:ly+1
             x0=[x0; (i-1)*lengthx/lx      -(j-1)*lengthy/ly];
         end
    else
        for j=1:(ly/2+1)
            x0=[x0; (i-1)*lengthx/lx      -(j-1)*lengthy/ly];
        end
    end
end
nodes=[];
for i=1:lx
    if i<lx/2+1
        for j=1:ly
            nodes=[nodes; (ly+1)*(i-1)+j (ly+1)*i+j (ly+1)*(i-1)+j+1;];
            nodes=[nodes; (ly+1)*i+j (ly+1)*i+j+1 (ly+1)*(i-1)+j+1;];
        end
    elseif i==lx/2+1
            for j=1:ly/2
                nodes=[nodes; (ly+1)*(i-1)+j (ly/2+1)*(i-lx/2-1)+j+153  (ly+1)*(i-1)+j+1;];
                nodes=[nodes; (ly/2+1)*(i-lx/2-1)+j+153 (ly/2+1)*(i-lx/2-1)+j+154 (ly+1)*(i-1)+j+1;];
            end
    else
        for j=1:ly/2
            nodes=[nodes; (ly/2+1)*(i-lx/2-2)+j+153 (ly/2+1)*(i-lx/2-1)+j+153 (ly/2+1)*(i-lx/2-2)+j+154;];
            nodes=[nodes; (ly/2+1)*(i-lx/2-1)+j+153 (ly/2+1)*(i-lx/2-1)+j+154 (ly/2+1)*(i-lx/2-2)+j+154;];
        end
    end
end
bcdof=[];
bcval=[];
for i=1:ly+1
        bcdof=[bcdof 1+2*(i-1) 2+2*(i-1)];
        bcval=[bcval  0   0];
end
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
matmtx=fematiso(1,emodule,poisson);     %constitutive matrice
for iel=1:nel       %loop for the total number of element
    nd(1)=nodes(iel,1);         %1st connected node for (iel)-th element
    nd(2)=nodes(iel,2);         %2nd connected node for (iel)-th element
    nd(3)=nodes(iel,3);         %3rd connected node for (iel)-th element
    x1=x0(nd(1),1); y1=x0(nd(1),2);     %coord values of 1st node
    x2=x0(nd(2),1); y2=x0(nd(2),2);     %coord values of 2nd node
    x3=x0(nd(3),1); y3=x0(nd(3),2);     %coord values of 3rd node
    index=feeldof(nd,nnel,ndof);        %extract system dofs for the element    
    area=abs(0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2)); %area of triangula
	area2=area*2;
     dhdx=(1/area2)*[(y2-y3) (y3-y1) (y1-y2)];       %derivatives w.r.t x
     dhdy=(1/area2)*[(x3-x2) (x1-x3) (x2-x1)];       %derivatives w.r.t y
 	kinmtx2=fekine2d(nnel,dhdx,dhdy);               %kinematic matrice   
   	k=kinmtx2'*matmtx*kinmtx2*area*h;              %element stiffness matrice
   	kk=feasmb_2(kk,k,index);                               %assemble element matrics
end
kk1=kk;
%--------------------------------------------------------------------------
%   force vector
%--------------------------------------------------------------------------
ff(442,1)=fload;
ff(137,1)=2*fload;
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
Patch_xy = zeros(6,nel);
for ielp=1:nel           % loop for the total number of elements
     nd(1)=nodes(ielp,1); % 1st connected node for (iel)-th element
     nd(2)=nodes(ielp,2); % 2nd connected node for (iel)-th element
     nd(3)=nodes(ielp,3); % 3rd connected node for (iel)-th element

    x1=x0(nd(1),1); y1=x0(nd(1),2);% coord values of 1st node
    x2=x0(nd(2),1); y2=x0(nd(2),2);% coord values of 2nd node
    x3=x0(nd(3),1); y3=x0(nd(3),2);% coord values of 3rd node
    Patch_xy(:,ielp)=[x1;x2;x3;y1;y2;y3];
    xcentre=(x1+x2+x3)/3; ycentre=(y1+y2+y3)/3;
    
    index=feeldof(nd,nnel,ndof);% extract system dofs associated with element
    %-------------------------------------------------------
    %  extract element displacement vector
    %-------------------------------------------------------
   for i=1:edof
        eldisp(i)=disp(index(i));
    end
    area=abs(0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2));  % area of triangule
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
%-------------------------------------------------------------
% Data output
%-------------------------------------------------------------
fid_out=fopen('result_beam02.plt','w');

fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
fprintf(fid_out,'VARIABLES="x" "y" "u" "v" "sigax"  "sigmay" "sigmaxy"\n');
fprintf(fid_out,'ZONE T="flow-field", N= %8d,E=%8d,ET=TRIANGLE, F=FEPOINT\n',nnode,nel);
disp=full(disp);
for i=1:nnode
       fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',x0(i,1),x0(i,2),disp(2*i-1),disp(2*i),stress_node(1,i),stress_node(2,i),stress_node(3,i));
end
for i=1:nel
      fprintf(fid_out,'%8d%8d%8d\n',nodes(i,1),nodes(i,2),nodes(i,3));
end
scatter(x0(:,1), x0(:,2),50,disp(1:2:end),"filled");
colorbar
clf
hold on
patch(Patch_xy(1:3,:), Patch_xy(4:6,:),stress.');
colorbar