clear
clc

% 定义单元数、节点数、每个单元的节点数、每个节点的自由度
nel=25;
nnode=12;
nnel=2;
ndof=2;
sdof=nnode*ndof;

% 设定节点坐标
gcoord(1,1)=15.24;gcoord(1,2)=15.24;
gcoord(2,1)=30.48;gcoord(2,2)=15.24;
gcoord(3,1)=45.72;gcoord(3,2)=15.24;
gcoord(4,1)=60.96;gcoord(4,2)=15.24;
gcoord(5,1)=76.2;gcoord(5,2)=15.24;
gcoord(6,1)=91.44;gcoord(6,2)=0.0;
gcoord(7,1)=76.2;gcoord(7,2)=0.0;
gcoord(8,1)=60.96;gcoord(8,2)=0.0;
gcoord(9,1)=45.72;gcoord(9,2)=0.0;
gcoord(10,1)=30.48;gcoord(10,2)=0.0;
gcoord(11,1)=15.24;gcoord(11,2)=0.0;
gcoord(12,1)=0.0;gcoord(12,2)=0.0;

% 设定每个单元的两个节点
nodes(1,1)=1;nodes(1,2)=2;
nodes(2,1)=2;nodes(2,2)=3;
nodes(3,1)=3;nodes(3,2)=4;
nodes(4,1)=4;nodes(4,2)=5;
nodes(5,1)=12;nodes(5,2)=11;
nodes(6,1)=11;nodes(6,2)=10;
nodes(7,1)=10;nodes(7,2)=9;
nodes(8,1)=9;nodes(8,2)=8;
nodes(9,1)=8;nodes(9,2)=7;
nodes(10,1)=7;nodes(10,2)=6;
nodes(11,1)=1;nodes(11,2)=11;
nodes(12,1)=2;nodes(12,2)=10;
nodes(13,1)=3;nodes(13,2)=9;
nodes(14,1)=4;nodes(14,2)=8;
nodes(15,1)=5;nodes(15,2)=7;
nodes(16,1)=1;nodes(16,2)=12;
nodes(17,1)=1;nodes(17,2)=10;
nodes(18,1)=2;nodes(18,2)=11;
nodes(19,1)=2;nodes(19,2)=9;
nodes(20,1)=3;nodes(20,2)=10;
nodes(21,1)=3;nodes(21,2)=8;
nodes(22,1)=4;nodes(22,2)=9;
nodes(23,1)=4;nodes(23,2)=7;
nodes(24,1)=5;nodes(24,2)=8;
nodes(25,1)=5;nodes(25,2)=6;

% 设定材料参数
for i=1:4
    props(i,1)=2e11;
    props(i,2)=0.0004;
end
for i=5:10
    props(i,1)=2e11;
    props(i,2)=0.0008;
end
for i=11:15
    props(i,1)=2e11;
    props(i,2)=0.0006;
end
for i=16:25
    props(i,1)=2e11;
    props(i,2)=0.0005;
end

% 设定约束
bcdof(1)=12; bcval(1)=0;
bcdof(2)=16; bcval(2)=0;
bcdof(3)=20; bcval(3)=0;
bcdof(4)=23; bcval(4)=0;
bcdof(5)=24; bcval(5)=0;


ff=zeros(sdof,1);
kk=zeros(sdof,sdof);
index=zeros(nnel*ndof,1);
k=zeros(nnel*ndof,nnel*ndof);
eldisp=zeros(nnel*ndof,1);

% 设定受力
ff(1)=5000;
ff(18)=-5000;
ff(24)=0;

% 计算刚度矩阵
for iel=1:nel
    nd(1)=nodes(iel,1);
    nd(2)=nodes(iel,2);
    x1=gcoord(nd(1),1);
    y1=gcoord(nd(1),2);
    x2=gcoord(nd(2),1);
    y2=gcoord(nd(2),2);
    
    L=sqrt((x2-x1)^2+(y2-y1)^2);
    lij=(x2-x1)/L;
    mij=(y2-y1)/L;
    T=[lij mij 0 0;0 0 lij mij];
    E=props(iel,1);
    A=props(iel,2);
    ke=(E*A/L)*[1 -1;-1 1];
    k=T'*ke*T;
    
    [index]=feeldof(nd,nnel,ndof);
    kk=feasmb_k(kk,k,index);
end
[kk,ff]=feaply(kk,ff,bcdof,bcval);
KK=kk;

KK(:,24)=[];
KK(:,23)=[];
KK(:,20)=[];
KK(:,16)=[];
KK(:,12)=[];

KK(24,:)=[];
KK(23,:)=[];
KK(20,:)=[];
KK(16,:)=[];
KK(12,:)=[];

% 计算节点位移
disp=ff'/kk;

% 计算杆内应力   
for iel=1:nel
    nd(1)=nodes(iel,1);
    nd(2)=nodes(iel,2);  
    x1=gcoord(nd(1),1);
    y1=gcoord(nd(1),2);
    x2=gcoord(nd(2),1);
    y2=gcoord(nd(2),2);   
    L=sqrt((x2-x1)^2+(y2-y1)^2);
    lij=(x2-x1)/L;
    mij=(y2-y1)/L;
    T=[lij mij 0 0;0 0 lij mij;];   
    B=(1.0/L)*[-1 1]*T;
    for inode=1:nnel
        eldisp((inode-1)*ndof+1,1)=disp((nd(inode)-1)*ndof+1);
        eldisp((inode-1)*ndof+2,1)=disp((nd(inode)-1)*ndof+2);
    end
    estrain=B*eldisp;                      
    E=props(iel,1);
    estress=E*estrain;                      
    stress(iel)=estress(1);             
end