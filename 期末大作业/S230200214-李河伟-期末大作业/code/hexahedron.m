clear all;
format short;

dictionary='gongziliang.inp';%读取Inp文件

[gcoords,nodes,nnpe,ndpn,ElasticM,PoissonR,bc_dof,load_set,dimension,element_type,SecType,Thickness] = loadinp(dictionary);

 bcval=[];
%边界条件
c=length(bc_dof);
for i=1:c

bcval=[bcval 0];

end

nel=max(nodes(:,1));     %总单元数

%各方向高斯积分点
nglx=2;
ngly=2;
nglz=2;

[point2,weight2]=feglqd2(nglx,ngly,nglz);    %积分点和权重

nnel=nnpe;
ndof=ndpn;
nnode=length(gcoords);
sdof=nnode*ndof;
edof=nnel*ndof;
ff=sparse(sdof,1);			
kk=sparse(sdof,sdof);
k=sparse(edof,edof);
elastic=ElasticM;
poisson=PoissonR;
 matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson   0   0    0; 
   poisson  (1-poisson)   poisson   0   0    0;
   poisson  poisson  (1-poisson)    0   0    0;
   0    0    0    (1-2*poisson)/2   0    0;
   0    0    0    0    (1-2*poisson)/2   0;
   0    0    0    0    0   (1-2*poisson)/2];

x=[];
y=[];
z=[];

for i=1:nnode
    
     x=[x gcoords(i,2)]; 
     y=[y gcoords(i,3)]; 
     z=[z gcoords(i,4)]; 
    
end

for j=1:nel   %循环计算所有单元
   k=sparse(edof,edof);
    for i=1:nnel              %每个单元节点编号及坐标
        
       nd(i)=nodes(j,i+1);     %一次取一个单元的序号
       ndx(i)=x(1,nd(i))       %本次循环各点的坐标值，按照1到8排列
       ndy(i)=y(1,nd(i))
       ndz(i)=z(1,nd(i))
    end

     for intx=1:nglx
        s=point2(intx,1);           %sampling point in x-axis    选取积分点S值，并配置权数
        wtx=weight2(intx,1);        %weight in x-axis           
        for inty=1:ngly
            n=point2(inty,2);       %sampling point in y-axis   选取积分点N值，并配置权数
            wty=weight2(inty,2);    %weight in y-axis
            for intz=1:nglz
                c=point2(intz,3);
                wtz=weight2(intz,3);
                
                 [shapeq,dhdsq,dhdnq,dhdzq]=feisoq8(s,n,c); %六面体形函数 
                 
                 [jacob3]=fejacob3(nnel,shapeq,dhdsq,dhdnq,dhdzq,ndx,ndy,ndz);
                 
                 detjacob=det(jacob3);
                 
                 invjacob=inv(jacob3);
                 
                 [dhdx,dhdy,dhdz]=federiv3(nnel,dhdsq,dhdnq,dhdzq,invjacob);
                 kinmtx2=fekine3d(nnel,dhdx,dhdy,dhdz);  
                 k=k+kinmtx2'*matmtrx*kinmtx2*wtx*wty*wtz*detjacob;
            end
        end
     end
     
              index=feeldof(nd,nnel,ndof);        %extract system dofs for the element                    提取系统自由度
              kk=feasmb_2(kk,k,index);            %assemble element matrics                               组装刚度矩阵
    
    
end

singularcheck(kk)   
disp=sparse(sdof,1);
kk1=kk;      %储存刚度矩阵

o=length(load_set(:,1));
for i=1:o

ff(load_set(i,1),1)=load_set(i,2);

end


[kk,ff]=feaplyc2(kk,ff,bc_dof,bcval);
[LL UU]=lu(kk);
utemp=LL\ff;
disp=UU\utemp;

n=max([max(nodes(:,1)),max(nodes(:,2)),max(nodes(:,4)),max(nodes(:,3)),max(nodes(:,5)),max(nodes(:,6)),max(nodes(:,7)),max(nodes(:,8)),max(nodes(:,9))]);
fid = fopen('gongziliang.plt','w'); 
fprintf(fid,'TITLE="test case governed by poisson equation"\n');
fprintf(fid,'VARIABLES="X""Y""Z""U""V""W"\n');
fprintf(fid,'ZONE N=%8d,E=%8d,ET=BRICK,F=FEPOINT\n',n,nel);

for i = 1:n
    fprintf(fid,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',x(1,i),y(1,i),z(1,i),disp(3*i-2)+0,disp(3*i-1)+0,disp(3*i)+0);
end

for i=1:nel
    fprintf(fid,'%8d%8d%8d%8d%8d%8d%8d%8d\n',nodes(i,2),nodes(i,3),nodes(i,4),nodes(i,5),nodes(i,6),nodes(i,7),nodes(i,8),nodes(i,9));
end


