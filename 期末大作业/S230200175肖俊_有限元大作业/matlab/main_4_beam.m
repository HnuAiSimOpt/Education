%% 
% 本程序针对四边形等参单元求解梯形梁(受均布力)问题进行计算
%模型：一端受到固定约束，另一端自由的梯形梁。
%，计算位移、应力和应变这三个指标的分量。
% 后处理包括：
%            1.梯形梁初始图形figure（1）
%            2.变形后的变形图figure（2）
%            2.变形后坐标、应力、应变云图figure（3）~figure（10）
% 输出内容包括：
%            1.图形输出；
%            2.plt文件输出(位移、应力、应变、材料常量等)。
%% 
clear all;clc;clear;
format short
first_time=cputime; 

%梁的节点和单元参数
lengthx=12;             %长度方向长度
lengthy=6;              %宽度方向长度
lx=60;                  %长度方向单元数
ly=30;                  %宽度方向单元数
nel=lx*ly;            %总的单元数
nnel=4;                %每个单元的节点数
ndof=2;                 %每个节点的自由度
nnode=(lx+1)*(ly+1);    %系统总节点数    
sdof=nnode*ndof;        %系统总自由度
edof=nnel*ndof;        %每个单元的自由度

%材料参数
Emodule=1E6;                  %弹性模量E
Poisson=0.3;                %泊松比

%高斯积分点的选择（可分别输出1X1高斯积分、2X2高斯积分、3X3高斯积分的结果）
%此处选择3X3高斯积分
% nglx=1; ngly=1;       %(1X1高斯积分)
% nglx=2; ngly=2;       %(2X2高斯积分)
nglx=3; ngly=3;         %(3X3高斯积分)
nglxy=nglx*ngly;        %高斯积分点总数
%%
%节点坐标
x0=[];
for i=1:lx+1
    for j=1:ly+1
        x0=[x0; (i-1)*lengthx/lx  -0.5*lengthy*(1+(lx+1-i)/lx)*(1-(j+16)/ly)];
    end
end

%各单元的节点坐标
nodes=[];
for i=1:lx
    for j=1:ly
        nodes=[nodes; (ly+1)*(i-1)+j (ly+1)*i+j (ly+1)*i+j+1 (ly+1)*(i-1)+j+1;];
    end
end

%绘制初始位移图形
figure(1) 
axis off                           %取消对坐标轴的一切设置
axis equal                         %严格控制各坐标的分度使其相等
hold on
 for ie=1:nel                      %各单元
    for j=1:nnel+1                 %连接成环
        j1=mod(j-1,nnel)+1;
        xp(j)=x0(nodes(ie,j1),1);  %各单元节点x坐标
        yp(j)=x0(nodes(ie,j1),2);  %各单元节点y坐标
    end
       plot(xp,yp,'-b')   
 end
legend('变形前状态') 
%%
%进行矩阵的初始化
ff=sparse(sdof,1);			%力初始化
kk=zeros(sdof,sdof);	    %系统矩阵
U=sparse(sdof,1);		    %初始系统位移
u=sparse(edof,1);		    %初始单元位移
stress=zeros(nel,4,3);      %初始应力矩阵（selem表示单元数；4表示一个单元对应的是四个高斯积分点，3表示每一个高斯积分点对应的三个应力值）;
strain=zeros(nel,4,3);      %初始化应变矩阵
B=sparse(3,8);              %B矩阵初始化
D=sparse(3,3);			    %D矩阵初始化

%将1X1高斯积分、2X2高斯积分、3X3高斯积分的积分点值和权系数放在同一个矩阵中，方便选择提取
ss=[0 0 0;-1/sqrt(3) 1/sqrt(3) 0; -sqrt(0.6) 0 sqrt(0.6)];%高斯积分点值矩阵
tt=[2 0 0;1 1 0;5/9 8/9 5/9];                             %将1X1高斯积分、2X2高斯积分、3X3高斯积分权系数矩阵
intpoint=ss(nglx,1:nglx);
weight=tt(nglx,1:nglx);
D=Emodule/(1-Poisson^2)*[1 Poisson 0;Poisson 1 0;0 0 (1-Poisson)/2];

%进行单元的循环，完成单元刚度矩阵的计算；
for iel=1:nel                           
    for i=1:nnel
        nd(i)=nodes(iel,i);                %提取第iel个单元的第i个节点
        xcoord(i)=x0(nd(i),1);             %提取第i节点的x坐标值
        ycoord(i)=x0(nd(i),2);             %提取第i节点的y坐标值
    end
k=sparse(edof,edof);			           %单元矩阵的初始化
jacob=zeros(2,2);                          %雅克比矩阵的初始化

%采用高斯积分进行单元刚度矩阵的求解
 for intx=1:nglx
        x=intpoint(intx);                  %在x轴上的高斯点
        wx=weight(intx);                   %在x轴上的权系数
        for inty=1:ngly
            y=intpoint(inty);              %在y轴上的高斯点
            wy=weight(inty);               %在y轴上的权系数
            [shape,dr,ds]=dshape(x,y);     %计算形函数及其偏导
            jacob(1,1)=dr(1)*xcoord(1)+dr(2)*xcoord(2)+dr(3)*xcoord(3)+dr(4)*xcoord(4);
            jacob(1,2)=dr(1)*ycoord(1)+dr(2)*ycoord(2)+dr(3)*ycoord(3)+dr(4)*ycoord(4);
            jacob(2,1)=ds(1)*xcoord(1)+ds(2)*xcoord(2)+ds(3)*xcoord(3)+ds(4)*xcoord(4);
            jacob(2,2)=ds(1)*ycoord(1)+ds(2)*ycoord(2)+ds(3)*ycoord(3)+ds(4)*ycoord(4);
            detjacob=det(jacob);           %求解雅克比矩阵的行列式         
            invjacob=[jacob(2,2) -jacob(1,2);-jacob(2,1) jacob(1,1)]/detjacob; %雅克比矩阵求逆
            dd=invjacob*[dr(1) dr(2) dr(3) dr(4);ds(1) ds(2) ds(3) ds(4)];     %全局坐标系下的偏微分
            B=bjuzhen(nnel,dd);                                                %求解B矩阵
            k=k+B'*D*B*wx*wy*detjacob;                                         %单元刚度矩阵求解
            
        end
 end
   m=1;
     for i=1:nnel            %进行单元刚度矩阵在系统刚度矩阵中的变换
         T(m)=2*nd(i)-1;
         T(m+1)=2*nd(i);
          m=m+2;
     end    
      for i=1:edof           %进行刚度矩阵的变换
          for j=1:edof
              kk(T(i),T(j))=kk(T(i),T(j))+k(i,j);
         end
      end
end
%%
%力的加载
for i=1:lx+1
ff(2*(ly+1)*i,1)=-1000;        %施加1000N的均布载荷
end

%边界条件的处理
w=1:ly+1;
b=size(w,2);                   %提取受约束的节点编号
for i=1:b
    kk(2*w(i)-1,2*w(i)-1)=1e20;%置大数法（这样可以保证约束点处的位移为0）
    kk(2*w(i),2*w(i))=1e20;
end

%矩阵求解
U=kk\ff;

%单元应力计算
  stab=[ 1.8660   -0.5000   0.1340 -0.50000 ;
          -0.500    1.8660 -0.5000 0.1340 ;
          0.1340   -0.5000  1.8660 -0.5000 ;
          -0.5000    0.1340 -0.5000 1.8660];    %积分点处应力与节点处应力的转换系数
      
for iel=1:nel                                 
    for i=1:nnel
        nd(i)=nodes(iel,i);                     %提取第iel个单元的第i个节点
        xcoord(i)=x0(nd(i),1);                  %提取第i节点的x坐标值
        ycoord(i)=x0(nd(i),2);                  %提取第i节点的y坐标值
    end
    k=sparse(edof,edof);			       
    
%进行高斯积分
    t=0;
    gstress=zeros(3,4);                         %每一列代表一个高斯点上的应力
    gstrain=zeros(3,4);
    for intx=1:nglx
        x=intpoint(intx);                       %在x轴上的高斯点
        wx=weight(intx);                        %在x轴上的权系数
        for inty=1:ngly
            y=intpoint(inty);                   %在y轴上的高斯点
            wy=weight(inty);                    %在y轴上的权系数
            [shape,dr,ds]=dshape(x,y);          %计算形函数及其偏导
            jacob(1,1)=dr(1)*xcoord(1)+dr(2)*xcoord(2)+dr(3)*xcoord(3)+dr(4)*xcoord(4);
            jacob(1,2)=dr(1)*ycoord(1)+dr(2)*ycoord(2)+dr(3)*ycoord(3)+dr(4)*ycoord(4);
            jacob(2,1)=ds(1)*xcoord(1)+ds(2)*xcoord(2)+ds(3)*xcoord(3)+ds(4)*xcoord(4);
            jacob(1,2)=ds(1)*xcoord(1)+ds(2)*ycoord(2)+ds(3)*ycoord(3)+ds(4)*ycoord(4);
            detjacob=det(jacob);                %雅克比矩阵的行列式         
            invjacob=[jacob(2,2) -jacob(1,2);-jacob(2,1) jacob(1,1)]/detjacob; %雅克比矩阵求逆
            
            dd=invjacob*[dr(1) dr(2) dr(3) dr(4);ds(1) ds(2) ds(3) ds(4)];     %全局坐标系下的偏微分
                        
            B=bjuzhen(nnel,dd);                 %求解几何矩阵B
   
            u=[U(2*nd(1)-1);U(2*nd(1));U(2*nd(2)-1);U(2*nd(2));U(2*nd(3)-1);U(2*nd(3));U(2*nd(4)-1);U(2*nd(4))];%将每个单元各个节点位移集合           
            estrain=B*u;                        %计算单元应变
            estress=D*estrain;                  %计算单元应力          
            t=t+1;                              %t表示第t个积分点
            gstress(:,t)=estress;               %最终形成3X4矩阵，每一列表示一个高斯点的应力
            gstrain(:,t)=estrain;
        end
    end

%实现高斯点向节点应力的转换
   for i=1:3
        for j=1:4
           for n=1:4
                stress(iel,j,i)= stress(iel,j,i)+stab(j,n)*gstress(i,n);
                strain(iel,j,i)=strain(iel,j,i)+stab(j,n)*gstrain(i,n);
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

% 每个节点的应力值与应变值
stress_node=zeros(3,nnode);	
strain_node=zeros(3,nnode);
for inode=1:nnode
    numel= indneigh(inode);
    for i=1:numel
        ind_nel= neigh_node{inode}(i);
        ind_nod=neigh_node_ind{inode}(i);
        for j=1:3
            stress_node(j,inode)=stress_node(j,inode)+stress(ind_nel,ind_nod,j);
            strain_node(j,inode)=strain_node(j,inode)+strain(ind_nel,ind_nod,j);
        end
    end
    stress_node(:,inode)=stress_node(:,inode)/numel;
    strain_node(:,inode)=strain_node(:,inode)/numel;
end
%%
%输出
%输出节点坐标、位移、应力和应变信息
fid_out=fopen('result.plt','w');
fprintf(fid_out,'Title="test case governed by poisson equation"\n');
fprintf(fid_out,'lengthx=12;lengthy=6;lx=60;ly=30;Emodule=1E6;Poisson=0.3"\n');
fprintf(fid_out,'Variables="x" "y" "u" "v" "sigamx"  "sigmay" "sigmaxy" "epsilongx" "epsilongy""epsilongxy"\n');
fprintf(fid_out,'ZONE T="flow-field", N= %8d,E=%8d,ET=QUADRILATERAL, F=FEPOINT\n',nnode,nel);
for i=1:nnode
    fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',x0(i,1),x0(i,2),U(2*i-1),U(2*i),stress_node(1,i),stress_node(2,i),stress_node(3,i),strain_node(1,i),strain_node(2,i),strain_node(3,i));
end
for i=1:nel
    fprintf(fid_out,'%8d%8d%8d%8d\n',nodes(i,1),nodes(i,2),nodes(i,3),nodes(i,4));
end

%输出变形后位移图形
figure(2)
axis off%取消对坐标轴的一切设置
axis equal%严格控制各坐标的分度使其相等
hold on
 for ie=1:nel
    for j=1:nnel+1
        j1=mod(j-1,nnel)+1;
        xp(j)=x0(nodes(ie,j1),1);
        yp(j)=x0(nodes(ie,j1),2);
        xm(j)=xp(j)+U(2*nodes(ie,j1)-1);
        ym(j)=yp(j)+U(2*nodes(ie,j1));
    end
       plot(xm,ym,'-r')   
 end
legend('变形后状态')     %绘制出变形前后的图形

% 输出云图
figure(3)
plotu(nel,U,nodes,x0,nnode,1)            %输出u方向位移云图
figure(4)
plotu(nel,U,nodes,x0,nnode,2)            %输出v方向位移云图
figure(5)
plotstress(nel,stress_node,nodes,x0,1)   %输出x方向sigamx云图
figure(6)
plotstress(nel,stress_node,nodes,x0,2)   %输出y方向sigmay云图
figure(7)
plotstress(nel,stress_node,nodes,x0,3)   %输出sigmaxy云图
figure(8)
plotstrain(nel,strain_node,nodes,x0,1)   %输出x方向epsilongx云图
figure(9)
plotstrain(nel,strain_node,nodes,x0,2)   %输出y方向epsilongy云图
figure(10)
plotstrain(nel,strain_node,nodes,x0,3)   %输出剪切epsilongxy云图




