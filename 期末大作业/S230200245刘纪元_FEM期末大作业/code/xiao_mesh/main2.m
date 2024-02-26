%-------------------------------------------------------------------------
%Problem: Cantilever beam : quadrilateral element 
% Using FEM-Q4 element
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
format short e             % format是指数据类型
first_time=cputime; 
%---------------------------------
%input data for control parameters

nnel=4;                    %number of nodes per element
ndof=2;                    %number of dofs per node           
emodule=2.1*power(10,11);  %elastic modulus
poisson=0.3;               %Poisson's ratio
nglx=2; ngly=2;            %2x2 Gauss-Legendre quadrature
nglxy=nglx*ngly;           %number of sampling points per element
%--------------------------------------------
%导入网格节点信息
x0 = load('C:\Users\alein\Desktop\S200200239_FEM\code\xiao_mesh\txt\x0.txt');
%导入网格节点与单元的连接信息
nodes = load('C:\Users\alein\Desktop\S200200239_FEM\code\xiao_mesh\txt\nodes.txt');
[nnode w]=size(x0);        %nnode：get total number of nodes in system  
[nel q]=size(nodes);       %nel：get number of element
sdof=nnode*ndof;           %total system dofs
edof=nnel*ndof;            %degrees of freedom per element
x0 = x0(:,[2 3]);
nodes = nodes(:,2:5);
%**************************************************************************
% Output the mesh
iplot=1;                   %施加平面力载荷
if iplot==1
   figure(99)
   hold on
   axis off
   axis equal
%搜索每个单元四个节点对应的坐标
   for ie=1:nel
        for j=1:nnel+1
            j1=mod(j-1,nnel)+1;
            xp(j)=x0(nodes(ie,j1),1);
            yp(j)=x0(nodes(ie,j1),2);
        end
        plot(xp,yp,'-')
   end
%pause           
end
%----------------------------------
%initialization of matrices and vectors
%--------------------------------------

ff=sparse(sdof,1);          %system force vector
fff=full(ff);               %复原力向量矩阵
k=sparse(edof,edof);		%initialization of element matrix
kk=zeros(sdof,sdof);		%system matrix
disp=sparse(sdof,1);		%system displacement vector
eldisp=sparse(edof,1);		%element displacement vector
stress=zeros(nel,4,3);		%matrix containing stress components
strain=sparse(nel,4,3);		%matrix containing strain components
index=sparse(edof,1);		%index vector
kinmtx=sparse(3,edof);		%kinematic matrix(这里的kinmtx矩阵是单元刚度矩阵)
matmtx=sparse(3,3);			%constitutive matrix（本构矩阵：应力应变关系，即有限元中的D矩阵）

%-------------------------------------------------
%compute element matrices and vectors and assemble
%-------------------------------------------------

[point2,weight2]=feglqd2(nglx,ngly);    %sampling points & weights
matmtx=fematiso(1,emodule,poisson);     %constitutive matrice
for iel=1:nel                           %loop for the total number of element
    for i=1:nnel
        nd(i)=nodes(iel,i);             %extract nodes for (iel)-th element
        xcoord(i)=x0(nd(i),1);          %extract x value of the node
        ycoord(i)=x0(nd(i),2);          %extract y value of the node
    end
    [angleA,angleB,angleC,angleD]=cal_angle(xcoord,ycoord);
    if angleA >=pi | angleA >=pi | angleA >=pi | angleA >=pi 
        temp='quadrilateral element is concave - FEM can not solve'
        break;
    end
    k=sparse(edof,edof);			    %initialization of element matrix
    
    %---------------------
    %numerical integration  数值积分
    %---------------------
    for intx=1:nglx
        x=point2(intx,1);               %sampling point in x-axis
        wtx=weight2(intx,1);            %weight in x-axis
        for inty=1:ngly
            y=point2(inty,2);           %sampling point in y-axis
            wty=weight2(inty,2);        %weight in y-axis
            
            [shape,dhdr,dhds]=feisoq4(x,y);                 %compute shape functions and derivatives at sampling point
            
            jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  %compute Jacobian
            
            detjacob=det(jacob2);                           %determinant of Jacobian           
            
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
%导入集中载荷所施加的节点编号信息
f_n = load('C:\Users\alein\Desktop\S200200239_FEM\code\xiao_mesh\txt\fload_number.txt');
%在局部坐标系下，对该节点的y方向施加集中载荷
ff(2*f_n,1)=-100*power(10,3);
%------------------------
%input data for boundary conditions
%----------------------------------
bcdof=[];
bcval=[];        
bc1 = load('C:\Users\alein\Desktop\S200200239_FEM\code\xiao_mesh\txt\bc1_number.txt');
%--------------------------------------
%apply boundary condition
%------------------------
[m n]=size(bc1);
for i=1:m
    a=bc1(i);
    bcdof=[bcdof 2*a-1 2*a];            %给边界约束条件节点处的自由度编号
    bcval=[bcval 0 0];                  %给定相对应自由度位移变化量
end
[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);
%-------------------------
%solve the matrix equation
%-------------------------

disp=kk\ff;                             %求解位移矩阵
%[LL UU]=lu(kk);
%utemp=LL\ff;
%disp=UU\utemp;

solve_time = cputime-first_time;        %计算实际运算所花时间

EU=0.5*disp'*kk1*disp;
U_fem=disp;                             %位移矩阵
%-------------------------------------------------
%  element stress computation
%-------------------------------------------------
energy=0;
stab=[1.866 -0.5 0.134 -0.5;
      -0.5 1.866 -0.5 0.134;
      0.134 -0.5 1.866 -0.5;
      -0.5 0.134 -0.5 1.866];
  
%下面搜寻节点坐标的代码重复了，可以不用重复循环，直接截取前面所得到的节点坐标矩阵即可
for iel=1:nel                           %loop for the total number of element
    for i=1:nnel
        nd(i)=nodes(iel,i);             %extract nodes for (iel)-th element
        xcoord(i)=x0(nd(i),1);          %extract x value of the node
        ycoord(i)=x0(nd(i),2);          %extract y value of the node
    end
    [angleA,angleB,angleC,angleD]=cal_angle(xcoord,ycoord);
    if angleA >=pi | angleA >=pi | angleA >=pi | angleA >=pi
        temp='quadrilateral element is concave - FEM can not solve'
        break; 
    end    
    k=sparse(edof,edof);			    %initialization of element matrix   
    
    %---------------------
    %numerical integration
    %---------------------
    pk=0;
    gstress=zeros(3,4);
    for intx=1:nglx
        x=point2(intx,1);               %sampling point in x-axis
        wtx=weight2(intx,1);            %weight in x-axis
        for inty=1:ngly
            y=point2(inty,2);           %sampling point in y-axis
            wty=weight2(inty,2);        %weight in y-axis
            
            [shape,dhdr,dhds]=feisoq4(x,y);                 %compute shape functions and derivatives at sampling point
            
            jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  %compute Jacobian
            
            detjacob=det(jacob2);                           %determinant of Jacobian           
            
            invjacob=inv(jacob2);                           %inverse of Jacobian matrix
            
            [dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob);  %derivatives w.r.t physical coordinate
            
            kinmtx2=fekine2d(nnel,dhdx,dhdy);               %kinematic matrice
            
            index=feeldof(nd,nnel,ndof);                    %extract system dofs for the element
            
            for i=1:edof
                eldisp(i,1)=disp(index(i));
            end            
            estrain=kinmtx2*eldisp;                         % compute strains
            estress=matmtx*estrain;                         % compute stresses           
            energy=energy+0.5*estrain'*matmtx*estrain*wtx*wty*detjacob;
            pk=pk+1;
            gstress(:,pk)=estress;
        end
    end
    for i=1:3
        for j=1:4
            for n=1:4
                stress(iel,j,i)= stress(iel,j,i)+stab(j,n)*gstress(i,n);
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
% Data output   将所有信息写入到result_beamQ4_1.plt文件
%-------------------------------------------------------------
fid_out=fopen('C:\Users\alein\Desktop\S200200239_FEM\code\xiao_mesh\result_beamQ4_2.plt','w');
fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
fprintf(fid_out,'VARIABLES="x" "y" "u" "v" "sigax"  "sigmay" "sigmaxy"\n');
fprintf(fid_out,'ZONE T="flow-field", N= %8d,E=%8d,ET=QUADRILATERAL, F=FEPOINT\n',nnode,nel);
for i=1:nnode
    fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',x0(i,1),x0(i,2),disp(2*i-1),disp(2*i),stress_node(1,i),stress_node(2,i),stress_node(3,i));
end
for i=1:nel
    fprintf(fid_out,'%8d%8d%8d%8d\n',nodes(i,1),nodes(i,2),nodes(i,3),nodes(i,4));
end