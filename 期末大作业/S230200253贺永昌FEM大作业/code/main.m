clc;
clear;
Material = [200*10^9;10^-3;7.85*10^3]; %单元材料参数（分别为杨氏弹性模量，横截面面积，密度）
DYS = ones(1,100); %单元数
DYJD = [1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 19 20 21 22 23 24 25 26 28 29 30 31 32 33 34 35 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 19 20 21 22 23 24 25 26 1  2  3  4  5  6  7  8;
      2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 20 21 22 23 24 25 26 27 29 30 31 32 33 34 35 36 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 20 21 22 23 24 25 26 27 29 30 31 32 33 34 35 36 29 30 31 32 33 34 35 36]; %单元节点编号
JDZB = [0 0 0;0 0 1;0 0 2;0 0 3;0 0 4;0 0 5;0 0 6;0 0 7;0 0 8;
               1 0 0;1 0 1;1 0 2;1 0 3;1 0 4;1 0 5;1 0 6;1 0 7;1 0 8;
               1 1 0;1 1 1;1 1 2;1 1 3;1 1 4;1 1 5;1 1 6;1 1 7;1 1 8;
               0 1 0;0 1 1;0 1 2;0 1 3;0 1 4;0 1 5;0 1 6;0 1 7;0 1 8]; %节点坐标
L = eLength(DYJD,JDZB);%每个单元的长度
ZDY = size(DYJD,2);%总单元数
ZJD = size(JDZB,1);%总节点数 
JDZYD = size(JDZB,2);%节点自由度
JDYS = zeros(JDZYD,ZJD);%节点约束(0为没有约束，1为有约束)
JDYS(1,1) = 1;JDYS(2,1) = 1;JDYS(3,1) = 1;
JDYS(1,10) = 1;JDYS(2,10) = 1;JDYS(3,10) = 1;
JDYS(1,19) = 1;JDYS(2,19) = 1;JDYS(3,19) = 1;
JDYS(1,28) = 1;JDYS(2,28) = 1;JDYS(3,28) = 1;%约束的节点及约束方向
P = zeros(JDZYD,ZJD);%节点载荷
P(3,27) = -1000;%输入节点载荷
IS = zeros(ZDY,JDZYD*2);%节点全局编号信息
K = zeros(ZJD*JDZYD);%总刚
stabK = zeros(ZJD*2);%总体几何矩阵(此程序只能计算平面稳定性问题）
L_d = zeros(2,ZDY);%局部单元节点位移
L_F = zeros(2,ZDY);%局部单元节点力 
G_F = zeros(JDZYD*2,ZDY);%全局单元节点力
PP = zeros(JDZYD,ZJD);%结构节点力(将不同单元的同一个节点的单元节点力相加，抵消相互作用力)
SG = zeros(1,ZDY);%单元应力（按单元存储）


%计算总刚
for i=1:ZDY
    E = Material(1,DYS(i));A = Material(2,DYS(i));
    %局部坐标系单刚
    L_Ke = E*A/L(i) * [1 -1;-1 1];
    TK = Transformation(DYJD(:,i),JDZB,L(i));
    %全局坐标系单刚
    G_Ke = TK' * L_Ke * TK;
    %节点全局编号
    for m = 1:JDZYD*2
        if m <= JDZYD
            IS(i,m) = (DYJD(1,i) - 1) * JDZYD + m;
        else
            IS(i,m) = (DYJD(2,i) - 1) * JDZYD + m - JDZYD;
        end
    end
    %集成总刚
    for m = 1:JDZYD*2
        for n = 1:JDZYD*2
            K(IS(i,m),IS(i,n)) = K(IS(i,m),IS(i,n)) + G_Ke(m,n);
        end
    end
    fprintf("第%d号单元局部单刚\n\n",i);disp(L_Ke);
    fprintf("第%d号单元转换阵\n\n",i);disp(TK);
    fprintf("第%d号单元全局单刚\n\n",i);disp(G_Ke);
end
fprintf("总刚度阵\n\n");disp(K);

%约束处理
for j = 1:ZJD
    for i = 1:JDZYD
        if JDYS(i,j) == 1
            IS_v = (j - 1) * JDZYD + i;
            K(IS_v,IS_v) = 10^30;%主对角元置大数
        end
    end
end

%计算全局节点位移
G_d = Cholesky(K,P,JDZYD);


%计算局部节点位移/局部节点力/单元应力/节点力
for i = 1:ZDY
    TK = Transformation(DYJD(:,i),JDZB,L(i));
    if JDZYD == 2
        E_d = [G_d(IS(i,1)) G_d(IS(i,2)) G_d(IS(i,3)) G_d(IS(i,4))];%i号单元全局节点位移（平面下）
    else 
        E_d = [G_d(IS(i,1)) G_d(IS(i,2)) G_d(IS(i,3)) G_d(IS(i,4)) G_d(IS(i,5)) G_d(IS(i,6))];%i号单元全局节点位移（空间下）
    end
    L_d(:,i) = TK * E_d';%i号局部单元节点位移
    E = Material(1,DYS(i));A = Material(2,DYS(i));
    L_Ke = E*A/L(i) * [1 -1;-1 1];
    L_F(:,i) = L_Ke * L_d(:,i);%i号局部单元节点力
    SG(i) = L_F(2,i) / A;%i号单元应力(第二个节点的局部节点力正负能反映拉正压负）
    G_F(:,i) = TK' * L_F(:,i);%i号单元全局节点力
    if JDZYD == 2
        s1 = [G_F(1,i) G_F(2,i)];s2 = [G_F(3,i) G_F(4,i)];
    else
        s1 = [G_F(1,i) G_F(2,i) G_F(3,i)];s2 = [G_F(4,i) G_F(5,i) G_F(6,i)];
    end
    PP(:,DYJD(1,i)) = PP(:,DYJD(1,i)) + s1';%将i号单元第一个节点的单元节点力存入PP列阵
    PP(:,DYJD(2,i)) = PP(:,DYJD(2,i)) + s2';%将i号单元第二个节点的单元节点力存入PP列阵
end

%计算支座反力
FR = PP - P;

%输出支座反力
fprintf("约束反力FR\n\n");
for j = 1:ZJD
    for i = 1:JDZYD
        if JDYS(i,j) == 1 %根据节点约束信息来判断那些节点施加了约束
            if i == 1
                fprintf("节点%dx方向约束反力",j);
                disp(FR(i,j));
            elseif i == 2
                fprintf("节点%dy方向约束反力",j);
                disp(FR(i,j));
            else 
                fprintf("节点%dz方向约束反力",j);
                disp(FR(i,j));
            end
        end
    end
end

%输出局部单元节点力
fprintf("\n\n局部单元节点力L_F\n\n");
disp(L_F);%矩阵的第二行正负代表实际杆件拉压状态

%输出局部单元节点位移
fprintf("\n\n局部单元节点位移L_d\n\n");
disp(L_d);%矩阵的第二行正负代表实际位移的正负

%输出节点位移
fprintf("\n\n全局单元节点位移G_d\n\n");
disp(G_d);

%输出单元内力
fprintf("\n\n单元内力\n\n");
disp(L_F(2,:));%单元内力其实就是局部单元节点力第二行的值

%输出单元应力
fprintf("\n\n单元应力SG\n\n");
disp(SG);


