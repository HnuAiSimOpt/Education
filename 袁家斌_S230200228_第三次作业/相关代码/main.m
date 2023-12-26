%姓名：袁家斌  学号：S230200228 
%purpose：用线性三角形单元对薄板结构进行应力应变分析
%边界条件：位移边界条件：节点1，4位移为0   力边界条件：节点2，3在y方向上力为0。
clc
clear;
numele=2;%单元数
numnode=4;%单元节点数
elenodecorre=[1 3 4;1 2 3];%各个单元分别对应的节点
%先将受载荷前的图画出来
ininodevector=[0 0 0.5 0 0.5 0.25 0 0.25];%各个节点的坐标
figure
axis off
axis equal
hold on
for inicoori=1:numnode
    plot(ininodevector(2*inicoori-1),ininodevector(2*inicoori),'bo');
end
for inicoorj=1:numele
    hold on
    line([ininodevector(2*elenodecorre(inicoorj,1)-1) ...
        ininodevector(2*elenodecorre(inicoorj,2)-1)],...
        [ininodevector(elenodecorre(inicoorj,1)*2) ...
        ininodevector(elenodecorre(inicoorj,2)*2)]);
    hold on
    line([ininodevector(2*elenodecorre(inicoorj,2)-1) ...
        ininodevector(2*elenodecorre(inicoorj,3)-1)],...
        [ininodevector(elenodecorre(inicoorj,2)*2) ...
        ininodevector(elenodecorre(inicoorj,3)*2)]);
    hold on
    line([ininodevector(2*elenodecorre(inicoorj,3)-1) ...
        ininodevector(2*elenodecorre(inicoorj,1)-1)],...
        [ininodevector(elenodecorre(inicoorj,3)*2) ...
        ininodevector(elenodecorre(inicoorj,1)*2)]);
end
%先给定相关参数
E=210e6;%弹性模量
NU=0.3;%泊松比
t=0.025;%薄板厚度；
k1=LinearTriangleElementStiffness(E,NU,t,0,0,0.5,0.25,0,0.25,1);%分别计算两个三角形单元的刚度矩阵k1和k2
k2=LinearTriangleElementStiffness(E,NU,t,0,0,0.5,0,0.5,0.25,1);
K=zeros(8,8);%因为有四个节点，整体的刚度矩阵为8x8，所以先生成一个8x8的空矩阵
K=LinearTriangleAssemble(K,k1,1,3,4);
K=LinearTriangleAssemble(K,k2,1,2,3);
%接下来要进行边界条件引入
%U1x=U1y=U4x=U4y=0----位移边界条件，1，4节点位移为0
%F2x=F3x=9.375;F2y=F3y=0----力边界条件；
k=K(3:6,3:6);
f=[9.375;0;9.375;0];
u=k\f;%用高斯消去法求解结点2，3位移
%得到节点2，3位移后求解节点1，4位移反力
U=[0;0;u;0;0];
F=K*U;
u1=[U(1);U(2);U(5);U(6);U(7);U(8)];
u2=[U(1);U(2);U(3);U(4);U(5);U(6)];
sigma1=LinearTriangleElementStresses(E,NU,t,0,0,0.5,0.25,0,0.25,1,u1);
sigma2=LinearTriangleElementStresses(E,NU,t,0,0,0.5,0,0.5,0.25,1,u2);%得到单元1和2的节点应力值
s1=LinearTriangleElementPStresses(sigma1);
s2=LinearTriangleElementPStresses(sigma2);
%相关的计算结果都在右边的对应表格中
finnodevector=[0 0 0.5+u(1,1) 0+u(2,1) 0.5+u(3,1) 0.25-u(4,1) 0 0.25];%受到载荷各个节点的坐标
hold on
%绘制受载荷后的图，由于变化太小，图中显示不明显
for fincoori=1:numnode
    plot(finnodevector(2*fincoori-1),finnodevector(2*fincoori),'ro');
end
for fincoorj=1:numele
    hold on
    line([finnodevector(2*elenodecorre(fincoorj,1)-1) ...
        finnodevector(2*elenodecorre(fincoorj,2)-1)],...
        [finnodevector(elenodecorre(fincoorj,1)*2) ...
        finnodevector(elenodecorre(fincoorj,2)*2)],'Color','red','LineStyle','--');
    hold on
    line([finnodevector(2*elenodecorre(fincoorj,2)-1) ...
        finnodevector(2*elenodecorre(fincoorj,3)-1)],...
        [finnodevector(elenodecorre(fincoorj,2)*2) ...
        finnodevector(elenodecorre(fincoorj,3)*2)],'Color','red','LineStyle','--');
    hold on
    line([finnodevector(2*elenodecorre(fincoorj,3)-1) ...
        finnodevector(2*elenodecorre(fincoorj,1)-1)],...
        [finnodevector(elenodecorre(fincoorj,3)*2) ...
        finnodevector(elenodecorre(fincoorj,1)*2)],'Color','red','LineStyle','--');
end