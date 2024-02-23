clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%导入数据
load('x0.mat');         %读取节点坐标x0为1225x2矩阵
load('nodes.mat');      %读取单元编号nodes为1152x4矩阵

nel=1152;               % 单元数
nnel=4;                 %每个单元节点数
ndof=2;                 %每个节点自由度
nnode=1225;             %节点总数   
sdof=2450;              %总自由度
edof=8;                 %每个单元自由度

%%%%%%%%%%%%%%%%%%%%%%%%材料属性

elastic=2e5;             %弹性模量
poisson=0.3;             %泊松比

fload=-10;               % 力10N
    
%%%%%%%%%%%%%%%%%%%%%%%%%%边界条件

bcdof=[];
bcval=[];
for i=1:25
        bcdof=[bcdof 1+2*(i-1) 2+2*(i-1)];
        bcval=[bcval 0 0];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%初始化矩阵
ff=sparse(sdof,1);			%力矩阵
k=sparse(edof,edof);		%单元刚度矩阵
kk=zeros(sdof,sdof);		%总刚度矩阵
disp=sparse(sdof,1);		%总结点位移
eldisp=zeros(edof,1);		%单元节点位移
stress=zeros(nel,4,3);		%总应力矩阵
strain=zeros(nel,4,3);		%总应变矩阵
B=sparse(3,8);		        % 应变矩阵B
D=sparse(3,3);			    % 材料矩阵D

D=elastic/(1-poisson*poisson)* ...
   [1  poisson 0; ...
   poisson  1  0; ...
   0  0  (1-poisson)/2];            %计算材料矩阵D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2x2高斯积分及总刚度矩阵组装
 
for iel=1:nel                       %遍历所有单元
   k=sparse(edof,edof); 
    for i=1:4
        nd(i)=nodes(iel,i);         %提取单元节点
        xcoord(i)=x0(nd(i),1);      %提取节点x坐标
        ycoord(i)=x0(nd(i),2);      %提取节点y坐标
    end
    pt=1/sqrt(3);
    gp = [-pt,-pt; pt, -pt; pt,pt; -pt,pt];             % 高斯积分点  
    w1= [1,1,1,1];
    for t = 1:4                                         % 遍历高斯积分点
            
        [shape,dhdr,dhds]=feisoq4(gp(t,1),gp(t,2));     %计算形函数以及对高斯计分点的偏导
            
         jacob=[dhdr;dhds]*[xcoord;ycoord]';            %计算雅克比矩阵  
            
         detjacob=det(jacob);                           %计算雅克比矩阵行列式        
            
         dxdy=(jacob)^(-1)*[dhdr;dhds];                 %计算形函数对节点偏导
         dhdx=dxdy(1,:);
         dhdy=dxdy(2,:);  
            
         B=[dhdx(1) 0 dhdx(2) 0 dhdx(3) 0 dhdx(4) 0;
                 0 dhdy(1) 0 dhdy(2) 0 dhdy(3) 0 dhdy(4);
                 dhdy(1) dhdx(1) dhdy(2) dhdx(2) dhdy(3) dhdx(3) dhdy(4) dhdx(4)];%计算矩阵B
                                                         
         k=k+B'*D*B*w1(t)*detjacob;                                               %计算单元刚度矩阵
            
    end
 
    G=zeros(8,2*nnode);                             %初始总刚变换矩阵
    G(1,2*nodes(iel,1)-1)=1;                        %总刚变换矩阵
    G(2,2*nodes(iel,1))=1;
    G(3,2*nodes(iel,2)-1)=1;
    G(4,2*nodes(iel,2))=1;
    G(5,2*nodes(iel,3)-1)=1;
    G(6,2*nodes(iel,3))=1;
    G(7,2*nodes(iel,4)-1)=1;
    G(8,2*nodes(iel,4))=1;
    kk=kk+G'*k*G;                                    %总体刚度矩阵合成
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%加载
for i=50:50:2450;
    ff(i,1)=fload;  
end                                                  %加均布载荷10N
%------------------------
%应用边界条件
[kk,ff]=feaplyc(kk,ff,bcdof,bcval);

[LL,UU]=lu(kk);                                     %杜立特分解
utemp=LL\ff;
disp=UU\utemp;                                      %节点位移

%********************计算节点应变应力*******************************

 stab=zeros(4);                                     %积分点应力到节点应力的转换矩阵
 for t=1:4
 stab(t,1)=0.25*(1-gp(t,1))*(1-gp(t,2));
 stab(t,2)=0.25*(1+gp(t,1))*(1-gp(t,2));
 stab(t,3)=0.25*(1+gp(t,1))*(1+gp(t,2));
 stab(t,4)=0.25*(1-gp(t,1))*(1+gp(t,2));
 end

 
 for iel=1:nel                      %%遍历所有单元
    for i=1:nnel
        nd(i)=nodes(iel,i);         %提取单元节点
        xcoord(i)=x0(nd(i),1);      %提取节点x坐标
        ycoord(i)=x0(nd(i),2);      %提取节点y坐标
    end
  
    k=sparse(edof,edof);			%初始单元刚度矩阵 
    
 
    %%%%%%%%%%%%%%2x2高斯积分
    pk=0;
    gstress=zeros(3,4);
    gstrain=zeros(3,4);
    pt=1/sqrt(3);
    gp = [-pt,-pt; pt, -pt; pt,pt; -pt,pt];             
    w1  = [1,1,1,1];
    for t = 1:length(w1)                                 % 遍历高斯积分点
            
        [shape,dhdr,dhds]=feisoq4(gp(t,1),gp(t,2));     %计算形函数以及对高斯计分点的偏导
            
         jacob=[dhdr;dhds]*[xcoord;ycoord]';           %计算雅克比矩阵 
            
         detjacob=det(jacob);                           %计算雅克比矩阵行列式        
            
         dxdy=(jacob)^(-1)*[dhdr;dhds];                 %计算形函数对节点的偏导
         dhdx=dxdy(1,:);
         dhdy=dxdy(2,:);  
            
         B=[dhdx(1) 0 dhdx(2) 0 dhdx(3) 0 dhdx(4) 0;
                  0 dhdy(1) 0 dhdy(2) 0 dhdy(3) 0 dhdy(4);
                  dhdy(1) dhdx(1) dhdy(2) dhdx(2) dhdy(3) dhdx(3) dhdy(4) dhdx(4)];%计算矩阵B
            
         eldisp=[disp(2*nodes(iel,1)-1);disp(2*nodes(iel,1));
                 disp(2*nodes(iel,2)-1);disp(2*nodes(iel,2));
                 disp(2*nodes(iel,3)-1);disp(2*nodes(iel,3));
                 disp(2*nodes(iel,4)-1);disp(2*nodes(iel,4))]; % 提取单元位移         
         estrain=B*eldisp;                                     % 计算高斯积分点的应变
         estress=D*estrain;                                    % 计算高斯积分点的应力         
         pk=pk+1;
         gstress(:,pk)=estress;                                 % 储存高斯积分点的应力
         gstrain(:,pk)=estrain;
    end
    for i=1:3                                                   % 节点三个方向应力
        for j=1:4                                               % 单元的四个节点应力
             stress(iel,j,i)= stress(iel,j,i)+stab(j,1)*gstress(i,1)+stab(j,2)*gstress(i,2)+stab(j,3)*gstress(i,3)+stab(j,4)*gstress(i,4);
             strain(iel,j,i)= strain(iel,j,i)+stab(j,1)*gstrain(i,1)+stab(j,2)*gstrain(i,2)+stab(j,3)*gstrain(i,3)+stab(j,4)*gstrain(i,4);
        end
    end
    
 end
stress_node=tolstress(nel,nnode,nodes,stress); %把stress转化为stress_node3x1225矩阵，3个方向应力1225个节点
strain_node=tolstress(nel,nnode,nodes,strain); %把strain转化为strain_node3x1225矩阵，3个方向应变1225个节点

%********************有限元后处理*******************************
% 输出网格和变形图
   figure(1)
   hold on
   axis off
   axis equal
   for ie=1:nel             
        for j=1:5
            j1=mod(j-1,4)+1;  
            xp(j)=x0(nodes(ie,j1),1);
            xp1(j)=x0(nodes(ie,j1),1)+disp((2*nodes(ie,j1)-1));
            yp(j)=x0(nodes(ie,j1),2);
            yp1(j)=x0(nodes(ie,j1),2)+disp(2*nodes(ie,j1));
        end
        plot(xp,yp,'-')
        hold on
        plot(xp1,yp1,'-g')
       
   end
% 输出云图
plotdisp(nel,disp,nodes,x0,1)            %输出u方向位移云图
plotdisp(nel,disp,nodes,x0,2)            %输出u方向位移云图
plotstrain(nel,strain_node,nodes,x0,1)   %输出x方向应变云图
plotstrain(nel,strain_node,nodes,x0,2)   %输出y方向应变云图
plotstrain(nel,strain_node,nodes,x0,3)   %输出剪切应变云图
plotstress(nel,stress_node,nodes,x0,1);  %输出x方向应力云图
plotstress(nel,stress_node,nodes,x0,2);  %输出y方向应力云图
plotstress(nel,stress_node,nodes,x0,3);  %输出剪切应力云图


