function [K]=Solve_K(poisson,element_solid,node,E)
eles_num=size(element_solid,1);   % 单元数
nodes_num=size(node,1);   % 节点数
total_dof=nodes_num*3;   % 总自由度数
nel=8;node_dof=3;
temp_element_solid=element_solid;temp_element_solid(:,1)=[];
temp_node=node;temp_node(:,1)=[];
element_dof=nel*node_dof;   % 一个单元自由度数
matmtx=E/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson   0   0    0; 
   poisson  (1-poisson)   poisson   0   0    0;
   poisson  poisson  (1-poisson)    0   0    0;
   0    0    0    (1-2*poisson)/2   0    0;
   0    0    0    0    (1-2*poisson)/2   0;
   0    0    0    0    0   (1-2*poisson)/2];   % 材料矩阵
K=sparse(total_dof,total_dof);		% 初始化整体刚度矩阵
nx=2;
ny=2;
nz=2;
% 调取权重、高斯积分点坐标
[point,weight]=Feglqd(nx,ny,nz);
for i=1:eles_num      
    for j=1:nel    
        nod(j)=temp_element_solid(i,j);    % 节点编号
        x(j)=temp_node(nod(j),1);   % 节点x坐标
        y(j)=temp_node(nod(j),2);   % 节点y坐标
        z(j)=temp_node(nod(j),3);   % 节点z坐标
    end 
    Ke=sparse(element_dof,element_dof);	   % 初始化单元刚度矩阵  
 for ix=1:nx
        sx=point(ix,1);        % x方向高斯积分点   
        sw=weight(ix,1);       % x方向权重           
        for iy=1:ny
            mx=point(iy,2);      % y方向高斯积分点  
            mw=weight(iy,2);     % y方向权重
            for iz=1:nz
                tx=point(iz,3);     % z方向高斯积分点
                tw=weight(iz,3);    % z方向高斯积分点
        [J,B]=Jacob(sx,mx,tx,nel,x,y,z);    % 求解雅可比矩阵及B矩阵
        det_J=det(J);
        Ke=Ke+B'*matmtx*B*sw*mw*tw*det_J; 
            end
        end
  end
        index=Element_Dof(nod,nel,node_dof);    % 提取单元自由度序号 
        K(index,index)=K(index,index) + Ke;    % 组装单元刚度矩阵  
end