%有限元大作业，用六面体单元解决问题
%学号:S230200217 姓名:甘霏
clear
close all
clc
%---------------------
%  导入节点和单元信息
%---------------------
node=dlmread('node.txt',',');%读入节点信息
elements=dlmread('elements.txt',',');%读入单元信息
%---------------
%  输入控制参数
%---------------
x0=node(:,2:4);%节点坐标
nodes=elements(:,2:9);%单元编码
nnode=size(x0,1);%节点总数
nel=size(nodes,1);%单元总数
nnel=size(nodes,2);%每个单元节点数
ndof=3;%每个节点自由度数
sdof=nnode*ndof;%系统总自由度数
edof=nnel*ndof;%每个单元自由度数
elastic=2.1e11;%杨氏模量
poisson=0.3;%泊松比
ng=2;%高斯积分点数
%-------------------
%  初始化矩阵和向量
%-------------------
F=zeros(sdof,1);%初始化节点力
kk=sparse(sdof,sdof);%总体刚度矩阵的空矩阵
stress=zeros(nnode,7);%初始化应力矩阵
ninode=zeros(nnode,1);%初始化节点次数矩阵
mistress=zeros(nnel,1);%初始化Mises应力矩阵
X=zeros(nnel,1);%初始化X坐标矩阵
Y=zeros(nnel,1);%初始化Y坐标矩阵
Z=zeros(nnel,1);%初始化Z坐标矩阵
%-----------------
%  计算矩阵和向量
%-----------------
D=Dc(4,elastic,poisson);%弹性矩阵
[p,w]=gaussquadrature(ng);%高斯积分点坐标值和加权系数
F(3*1160-1)=1000;%1160点受y方向力
bcdof=[];%建立约束矩阵的空矩阵
for ibcdof=1:nnode
    if x0(ibcdof,3)==0
        bcdof=[bcdof ibcdof];%约束矩阵
    end
end
nbc=size(bcdof,2);%约束点数
for ini=1:nel
    for inj=1:nnel
        ninode(nodes(ini,inj))=ninode(nodes(ini,inj))+1;%每个节点用的次数
    end
end
%---------------------
%  计算刚度矩阵并组装
%---------------------
for iel=1:nel
    nd=nodes(iel,:);%第iel个单元包含的各节点
    for i=1:nnel
        X(i)=x0(nd(i),1);
        Y(i)=x0(nd(i),2);
        Z(i)=x0(nd(i),3);%第i个单元的各节点坐标
    end
    k=zeros(edof,edof);%初始化单元刚度矩阵
    %--------------------------
    %  用高斯积分求单元刚度矩阵
    %--------------------------
    for iz=1:ng
        z=p(iz);
        wz=w(iz);
        for iy=1:ng
            y=p(iy);
            wy=w(iy);
            for ix=1:ng
                x=p(ix);
                wx=w(ix);
                [N,dndr,dnds,dndt]=shape(x,y,z);%形函数
                J=Jacobain(nnel,dndr,dnds,dndt,X,Y,Z);%雅克比矩阵
                B=Bc(J,nnel,dndr,dnds,dndt);%单元应变矩阵
                k=k+(B')*D*B*wx*wy*wz*det(J);%计算单元刚度矩阵
            end
        end
    end
    %---------------
    %  组装刚度矩阵
    %---------------
    kk=zk(kk,nd,k,nnel,ndof);%总体刚度矩阵
end
%---------------
%  施加边界条件
%---------------
kk=boundary(bcdof,nbc,kk,ndof);
%---------------
%  计算节点位移
%---------------
disp=kk\F;
%---------------------
%  计算应力矩阵并组装
%---------------------
for iel=1:nel
    nd=nodes(iel,:);%第iel个单元包含的各节点
    for i=1:nnel
        X(i)=x0(nd(i),1);
        Y(i)=x0(nd(i),2);
        Z(i)=x0(nd(i),3);%第i个单元的各节点坐标
    end
    gstress=[];%初始化高斯点应力矩阵
    %----------------------------------
    %  用高斯积分点求高斯积分点应力矩阵
    %----------------------------------
    for iz=1:ng
        z=p(iz);
        for iy=1:ng
            y=p(iy);
            for ix=1:ng
                x=p(ix);
                [N,dndr,dnds,dndt]=shape(x,y,z);%形函数
                J=Jacobain(nnel,dndr,dnds,dndt,X,Y,Z);%雅克比矩阵
                B=Bc(J,nnel,dndr,dnds,dndt);%单元应变矩阵
                stress_g=ssc(disp,ndof,nd,B,D);%高斯点应力
                gstress=[gstress stress_g];%高斯点应力矩阵
            end
        end
    end
    %---------------------
    %  重组高斯积分点应力
    %---------------------
    gstress=[gstress(:,5),gstress(:,6),gstress(:,8),gstress(:,7),...
        gstress(:,1),gstress(:,2),gstress(:,4),gstress(:,3)];
    %-----------
    %  应力恢复
    %-----------
    a=sqrt(3);
    h=[recovery(-a,-a, a)
       recovery( a,-a, a)
       recovery( a, a, a)
       recovery(-a, a, a)
       recovery(-a,-a,-a)
       recovery( a,-a,-a)
       recovery( a, a,-a)
       recovery(-a, a,-a)];%应力恢复矩阵
    estress=h*transpose(gstress);%单元节点应力
    %---------------------------------
    %  计算Mises应力并重组单元节点应力
    %---------------------------------
    for mi=1:nnel
        s1=estress(mi,1);
        s2=estress(mi,2);
        s3=estress(mi,3);
        mistress(mi)=(((s1-s2)^2+(s2-s3)^2+(s3-s1)^2)/2)^0.5;%Mises应力
    end
    estress=[estress,mistress];%重组单元节点应力
    %-----------------
    %  组装节点总应力
    %-----------------
    for is=1:nnel
        stress(nd(is),:)=stress(nd(is),:)+estress(is,:);%所有单元叠加的节点应力
    end
end
%------------
%  应力平滑
%------------
for inode=1:nnode
    an=ninode(inode);
    stress(inode,:)=stress(inode,:)/an;%每个节点应力平滑
end
%--------------------------
%  将结果导出为文件jieguo.plt中
%--------------------------
nel = size(nodes,1);
node= size(x0,1);
nodedata=fopen('jieguo.plt','w');  
fprintf(nodedata,'TITLE="data"\n');
fprintf(nodedata,'VARIABLES=,"X", "Y","Z" ,"sig1","sig2","sig3","sig4","sig5","sig6","sig7"\n');
fprintf(nodedata,'ZONE T="%d  "  ,  N=%d, E=%d, ET=BRICK, F=FEPOINT\n',1,node,nel);
sigm1=stress(:,1);
sigm2=stress(:,2);
sigm3=stress(:,3);
sigm4=stress(:,4);
sigm5=stress(:,5);
sigm6=stress(:,6);
sigm7=stress(:,7);
for i=1:node
    fprintf(nodedata,'%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f\n',...
      x0(i,1),x0(i,2),x0(i,3),sigm1(i),sigm2(i),sigm3(i),sigm4(i),sigm5(i),sigm6(i),sigm7(i));
end  
for i=1:nel
    for j=1:size(nodes,2)
        fprintf(nodedata,'%d       ',nodes(i,j));  
    end
    fprintf(nodedata,'\n');
end