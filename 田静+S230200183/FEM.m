clear
clc
% 学号：S230200183      姓名：田静
% --------------------------------------------------------------------
% 就PPT里面例题，编写一个自己的解题代码。
% 采用两个三角形单元离散，单元1由节点1、2、3组成；单元2由节点1、3、4组成
% 考虑一个平面应力问题，如PPT案例所示，假设厚度h=1，材料为各项同性，杨氏模量为E=1，泊松比为ν=0，
% 自己代码计算结果
% 位移[0,0,0,0,-0.628175519630482,-7.68591224018472,1.25635103926096,-8.16628175519627]
% 应变：单元1[-0.157043879907621;0;-1.92147806004618]
%       单元2[0.314087759815241;-0.480369515011546;-0.157043879907620]
% 应力：单元1[-0.157043879907621;0;-0.960739030023090]
%       单元2[0.314087759815241;-0.480369515011546;-0.0785219399538102]
% 计算结果与PPT案例结果一致
% --------------------------------------------------------------------
% PPT中所给代码，该代码更加通用化，通过等分梯形x方向以及y方向长度，从而调整单元数量
% 令lx=1;ly=1;使得单元数量为2，计算位移结果与自己编写代码一致
% 此外，该代码在后处理计算应变、应力时，不仅仅给出了各个单元应变、应力，同时采用节点应力平均，获得节点应力
% element       1 :first_node:       1 second_node:       3 third_node:       2
% element       2 :first_node:       3 second_node:       4 third_node:       2
% node       1 :x:    0.000000e+00 y:   -2.000000e+00 u:    0.000000e+00 v:    0.000000e+00 sigmax:   -1.570439e-01 sigmay:    0.000000e+00 sigmaxy:   -9.607390e-01
% node       2 :x:    0.000000e+00 y:    0.000000e+00 u:    0.000000e+00 v:    0.000000e+00 sigmax:    7.852194e-02 sigmay:   -2.401848e-01 sigmaxy:   -5.196305e-01
% node       3 :x:    4.000000e+00 y:   -1.000000e+00 u:   -6.281755e-01 v:   -7.685912e+00 sigmax:    7.852194e-02 sigmay:   -2.401848e-01 sigmaxy:   -5.196305e-01
% node       4 :x:    4.000000e+00 y:    0.000000e+00 u:    1.256351e+00 v:   -8.166282e+00 sigmax:    3.140878e-01 sigmay:   -4.803695e-01 sigmaxy:   -7.852194e-02
% --------------------------------------------------------------------

%
% 自己编的代码
% -----------------------------------------------------
% 相关力和位移边界条件如图下所示：
% 载荷为节点4的集中载荷(方向向下)，位移约束为节点1、2处的固定约束
% 求问题各节点位移u、v和应力σx，σy和σxy
% 节点坐标如下：
% 节点1(0,0);节点2(0,-2);节点3(4,-1);节点4(4,0);
Element_number=2; %单元数量
No_nel=3;%每个单元节点数量
No_dof=2;%每个节点自由度
Node_number=4;%系统节点数量
Prop(1)=1;%弹性模量
Prop(2)=0;%泊松比
h=1;%厚度
%赋予节点坐标值
gcoord(1,1)=0;
gcoord(1,2)=0;
gcoord(2,1)=0;
gcoord(2,2)=-2;
gcoord(3,1)=4;
gcoord(3,2)=-1;
gcoord(4,1)=4;
gcoord(4,2)=0;
%赋予每个单元所包含节点编号
for No_element=1:2
     if No_element==1
          nodes(No_element,1)=1;
          nodes(No_element,2)=2;
          nodes(No_element,3)=3;
     else
          nodes(No_element,1)=1;
          nodes(No_element,2)=3;
          nodes(No_element,3)=4;
     end
end
%系统位移约束
ed(1:Node_number,1:2)=1;
constraint=[1 1;1 2;2 1;2 2];
for loopi=1:length(constraint)
    ed(constraint(loopi,1),constraint(loopi,2))=0;
end
%定义系统节点位移应属于哪一个自由度
dof=0;
for loopi=1:Node_number
    for loopj=1:2
        if ed(loopi,loopj)~=0
            dof=dof+1;
            ed(loopi,loopj)=dof;
        end
    end
end
%矩阵初始化
k=zeros(dof,dof);
f=zeros(Node_number*No_dof,1);
disp=zeros(dof,1);
eldisp=zeros(No_nel*No_dof,1);
e2s(1:6)=0;
f(8)=-1;
%系统刚度矩阵的组装
for loopi=1:Element_number
    for zi=1:3 %单元各节点对应自由度
        e2s((zi-1)*2+1)=ed(nodes(loopi,zi),1);
        e2s((zi-1)*2+2)=ed(nodes(loopi,zi),2);
    end
    if loopi==1
        x1=gcoord(1,1);y1=gcoord(1,2);
        x2=gcoord(2,1);y2=gcoord(2,2);
        x3=gcoord(3,1);y3=gcoord(3,2);
    else
        x1=gcoord(1,1);y1=gcoord(1,2);
        x2=gcoord(3,1);y2=gcoord(3,2);
        x3=gcoord(4,1);y3=gcoord(4,2);
    end
    area=0.5*(x1*y2+x2*y3*x3*y1-x1*y3-x2*y1-x3*y2);
    dNdx=(1/(2*area))*[(y2-y3) (y3-y1) (y1-y2)];
    dNdy=(1/(2*area))*[(x3-x2) (x1-x3) (x2-x1)];
    for i=1:3
        i1=(i-1)*2+1;
        i2=i1+1;
        B(1,i1)=dNdx(i);
        B(2,i2)=dNdy(i);
        B(3,i1)=dNdy(i);
        B(3,i2)=dNdx(i);
    end
    L_N{loopi}=B;
    D=Prop(1)/(1-Prop(2)*Prop(2))*[1 Prop(2) 0;Prop(2) 1 0;0 0 (1-Prop(2))/2];
    ke=h*area*B'*D*B;
    for jx=1:6
        for jy=1:6
            if (e2s(jx)*e2s(jy)~=0)
                k(e2s(jx),e2s(jy))=k(e2s(jx),e2s(jy))+ke(jx,jy);
            end
        end
    end
end
%去除外力中已约束位移的自由度
Number_con=-1;
for loopi=1:length(constraint)
    Number_con=Number_con+1;
    f((constraint(loopi,1)-1)*No_dof+constraint(loopi,2)-Number_con)=[];
end
%求解位移
disp=k\f;
%组装全部(包含约束位移)位移
flag=1;
iter=-1;
disp1=disp;
while flag
    for loopi=1:length(constraint)
        iter=iter+1;
        for loopj=1:(dof+iter)
            if loopj<((constraint(loopi,1)-1)*No_dof+constraint(loopi,2))
                dispt(loopj)=disp1(loopj);
            else
                dispt(loopj+1)=disp1(loopj);
            end
        end
        dispt((constraint(loopi,1)-1)*No_dof+constraint(loopi,2))=0;
        disp1=dispt;
    end
    if length(disp1)>=Node_number*No_dof
        flag=0;
    end
end
%后处理，计算应变、应力
for loopi=1:Element_number
    if loopi==1
       strain(:,loopi)=L_N{loopi}*disp1(:,1:6)';
       stress(:,loopi)=D*strain(:,loopi);
    else
       strain(:,loopi)=L_N{loopi}*[disp1(:,1:2) disp1(:,5:6) disp1(:,7:8)]';
       stress(:,loopi)=D*strain(:,loopi);
    end
end
%}

%{
% PPT中的代码
%------------------------------------------------------------------------------------------
% 以下是PPT中所给代码，该代码更加通用化，通过等分梯形x方向以及y方向长度，从而调整单元数量
% 令lx=1;ly=1;使得单元数量为2，计算位移结果与自己编写代码一致
% 该代码此外在后处理计算应变、应力时，不仅仅给出了各个单元应变、应力，同时采用节点应力平均，获得节点应力
%------------------------------------------------------------------------------------------
clear all
clc
first_time=cputime; 
format long
%--------------------------------------------
%input data for control parameters
%--------------------------------------------
lengthx=4;              %length of x-axis side of problem
lengthy=2;              %length of y-axis side of problem
Prop(1)=1;              %弹性模量
Prop(2)=0;              %泊松比
h=1;                    %厚度 
fload=-1;               % the total load
lx=1;                  % number of element in x-axis
ly=1;                   % number of element in y-axis
nel=2*lx*ly;            % number of element
nnel=3;                 %number of nodes per element
ndof=2;                 %number of dofs per node
nnode=(lx+1)*(ly+1);    %total number of nodes in system    
sdof=nnode*ndof;        %total system dofs
edof=nnel*ndof;         %degrees of freedom per element
%-------------------------------------------------------
%编制各节点对应坐标
%-------------------------------------------------------
x0=[];
for i=1:lx+1
    for j=1:ly+1
          x0=[x0; (i-1)*lengthx/lx,-0.5*lengthy*(1+(lx+1-i)/lx)*(1-(j-1)/ly)];
    end
end
%----------------------------------------------------------------------
%编制各单元对应节点编号
%----------------------------------------------------------------------
nodes=[];
for i=1:lx
    for j=1:ly
        nodes=[nodes;(ly+1)*(i-1)+j,(ly+1)*i+j,(ly+1)*(i-1)+j+1;];
        nodes=[nodes;(ly+1)*i+j,(ly+1)*i+j+1,(ly+1)*(i-1)+j+1;];
    end
end
%----------------------------------
%边界位移约束
%----------------------------------
bcdof=[];
bcval=[];
for i=1:ly+1
        bcdof=[bcdof 1+2*(i-1) 2+2*(i-1)];
        bcval=[bcval  0   0];
end
%-------------------------------------------------
%矩阵初始化
%-------------------------------------------------
ff=sparse(sdof,1);			%系统力矩阵
k=sparse(edof,edof);		%单元刚度矩阵
kk=sparse(sdof,sdof);		%系统刚度矩阵
disp=sparse(sdof,1);		%系统位移矩阵
eldisp=sparse(edof,1);		%单元位移矩阵
stress=zeros(nel,3);		%系统应力矩阵
strain=zeros(nel,3);		%系统应变矩阵
index=sparse(edof,1);		%index vector
kinmtx=sparse(3,edof);		%形函数对坐标求偏导L_N
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
    index=feeldof(nd,nnel,ndof);        %提取每个单元各节点对应的系统自由度 
    %---------------------------------------
    %形函数求偏导
    %---------------------------------------
    area=0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2); %area of triangula
    area2=area*2;	
    dhdx=(1/area2)*[(y2-y3) (y3-y1) (y1-y2)];       %derivatives w.r.t x
    dhdy=(1/area2)*[(x3-x2) (x1-x3) (x2-x1)];       %derivatives w.r.t y
    kinmtx2=fekine2d(nnel,dhdx,dhdy);            %L_N 
    D=Prop(1)/(1-Prop(2)*Prop(2))*[1 Prop(2) 0;Prop(2) 1 0;0 0 (1-Prop(2))/2];
    k=kinmtx2'*D*kinmtx2*area*h;              %单元矩阵
    kk=feasmb_2(kk,k,index);                     %刚度矩阵组装
end
kk1=kk;
%--------------------------------------------------------------------------
%施加力边界条件
%--------------------------------------------------------------------------
ff(sdof,1)=fload;
%------------------------
%施加位移边界条件
%------------------------
[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);
%-------------------------
%求解位移
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
    estress=D*estrain;             % compute stresses
    for i=1:3
        strain(ielp,i)=estrain(i);          % store for each element
        stress(ielp,i)=estress(i);          % store for each element          
    end
    energy=energy+0.5*estrain'*D*estrain*area*h;
end
neigh_node=cell(nnode,1);%每个节点分别有几个单元共用
indneigh=zeros(1,nnode);%每个节点所共用单元的编号
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
% 数据输出
%-------------------------------------------------------------
fid_out=fopen('result_beam01.plt','w');
for i=1:nel
      fprintf(fid_out,'element%8d :first_node:%8d second_node:%8d third_node:%8d\n',(i),nodes(i,1),nodes(i,2),nodes(i,3));
end
for i=1:nnode
      fprintf(fid_out,'node%8d :x:%16.6e y:%16.6e u:%16.6e v:%16.6e sigmax:%16.6e sigmay:%16.6e sigmaxy:%16.6e\n',(i),x0(i,1)+0,x0(i,2)+0,disp(2*i-1)+0,disp(2*i)+0,stress_node(1,i)+0,stress_node(2,i)+0,stress_node(3,i)+0);
end
type result_beam01.plt

%子函数
function [index]=feeldof(nd,nnel,ndof)
%----------------------------------------------------------
%  Purpose:
%     Compute system dofs associated with each element 
%     Variable Description:
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
%}





    
    
    



