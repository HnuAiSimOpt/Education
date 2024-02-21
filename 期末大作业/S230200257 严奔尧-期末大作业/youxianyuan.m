clc
clear
format compact
format shortG
jd=input('请输入节点数：');
dy=input('请输入单元数：');
E=input('请输入杨氏模量E：');
I=input('请输入惯性矩I：');
L=input('请输入单元长度L：');
A=input('请输入单元截面积：');
FAI=input('请输入单元相对旋转角度：');
%输入对应关系时，小节点放前面[单元 节点1 节点2]
dy_jd=input('请输入单元与节点对应关系：');
%输入力与扭矩约束[值 作用节点 作用类型]（转矩为3 x方向为1 y方向为2）
lys=input('力与转矩约束矩阵：');
%输入结构约束[作用节点 作用类型]（转角为3 x方向为1 y方向为2）
wys=input('结构约束矩阵：');
%原始数据
% L=1;
% E=3*10^10;
% P=1000;
% A=0.05;
% dy=2;jd=3;LL=[L 2*L];I=20*A;
% dy_jd=[1 1 2;2 2 3];
% FAI=[pi/2 0];
% q=P/L;M=P*L/10;
% lys=[44/125*P 1 1;-12*P*L/125 1 3;81/125*P 2 1;-P 2 2;-67/750*P*L 2 3;-P 3 2;P*L/3 3 3];
% wys=[1 1;1 2;1 3;3 1;3 2;3 3];
%对力约束与位移约束式子分别进行编号处理
wys(:,3)=(wys(:,1)-1)*3+wys(:,2);
lys(:,4)=(lys(:,2)-1)*3+lys(:,3);
%对力约束与位移约束式子进行排序
lys=sortrows(lys,4);
wys=sortrows(wys,3);
%单元刚度矩阵
syms fai e a i l real
k=[e*a/l 0 0 -e*a/l 0 0;
    0 12*e*i/l^3 6*e*i/l^2 0 -12*e*i/l^3 6*e*i/l^2;
    0 6*e*i/l^2 4*e*i/l 0 -6*e*i/l^2 2*e*i/l;
    -e*a/l 0 0 e*a/l 0 0;
    0 -12*e*i/l^3 -6*e*i/l^2 0 12*e*i/l^3 -6*e*i/l^2;
    0 6*e*i/l^2 2*e*i/l 0 -6*e*i/l^2 4*e*i/l];
t=[  cos(fai), sin(fai), 0;    
    -sin(fai), cos(fai), 0;
    0,        0, 1];
%坐标变换矩阵
T=blkdiag(t,t);
%总体坐标系下的单元刚度矩阵
K=T'*k*T;
%带入每个单元的数，生成单元刚度矩阵kk，其每一页对应相应页数的单元的刚度矩阵
for j=1:dy;
    e=E;
    i=I;
    l=LL(j);
    a=A;
    fai=FAI(j);
    kk(:,:,j)=eval(K);
end
%生成总体刚度矩阵KK
%采用元胞数组的方式对各项进行保存
%生成空元胞数组，元胞数组的行列大小与节点数相同
for j=1:jd;
    for jj=1:jd;
        ling1{j,jj}=zeros(3);
    end
end
ling2=ling1;
%将对单元刚度矩阵部分分成4分加入元胞数组中
for j=1:dy;
    kk1=kk(1:3,1:3,j);
    kk2=kk(1:3,4:6,j);
    kk3=kk(4:6,1:3,j);
    kk4=kk(4:6,4:6,j);  
    ling2{dy_jd(j,2),dy_jd(j,2)}=kk1+ling2{dy_jd(j,2),dy_jd(j,2)};
    ling2{dy_jd(j,2),dy_jd(j,3)}=kk2+ling2{dy_jd(j,2),dy_jd(j,3)};
    ling2{dy_jd(j,3),dy_jd(j,2)}=kk3+ling2{dy_jd(j,3),dy_jd(j,2)};
    ling2{dy_jd(j,3),dy_jd(j,3)}=kk4+ling2{dy_jd(j,3),dy_jd(j,3)};
end
%将元胞数组进行拼接，形成总体刚度矩阵
for j=1:jd;
    ling3(:,:,j)=cat(2,ling2{j,:});
end
KK=ling3(:,:,1);
for j=2:jd;
    KK=[KK;ling3(:,:,j)];
end
%消去有已知位移的行与列
b=KK;
b(:,wys(:,3))=[];
b(wys(:,3),:)=[];
kjiejuzhen=inv(b);
%提取对应外力
lyss=lys;
for j=1:size(wys,1);
    for jj=1:size(lys,1);
        if lyss(jj,4)==wys(j,3);
           lyss(jj,:)=0;
        end
        if jj==size(lyss,1);
           break
        end
    end
end
lyss(all(lyss==0,2),:)=[];
%求解weiyijie=[作用值 作用节点 作用类型（转角为3 x方向为1 y方向为2） 序列]
weiyijie=kjiejuzhen*lyss(:,1);
weiyijie(:,1)=weiyijie;
weiyijie(:,2)=lyss(:,2);
weiyijie(:,3)=lyss(:,3);
weiyijie(:,4)=lyss(:,4);
%计算不计作用在约束方向上时的支反力lysjiee=[作用值 作用节点 作用类型（转角为3 x方向为1 y方向为2） 序列]
lysjie(:,1)=KK(wys(:,3),lyss(:,4))*weiyijie(:,1);
lysjie(:,2:4)=wys(:,1:3);
%将作用在约束方向上时的支反力加在上面的求解结果上
for j=1:size(lysjie,1)
    for jj=1:size(lys,1);
        if lysjie(j,4)==lys(jj,4);
            lysjie(j,1)=lysjie(j,1)-lys(jj,1);
        end
    end
end
%答案
weiyijie 
lysjie
