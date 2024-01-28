%将三角形薄板划分为三角形单元格网格，每个三角单元格高为2m
clear 
format short e 
EleNum = 36; %单元格个数
PointNum = 28; %结点个数
VFixNum = 2;%受约束边界点数
ForceNum = 1;%载荷数
Young = 2e11;%弹性模量
Poiss = 0.25;%泊松比
Thick = 0.2;%薄板厚度

ForcePos = [1 0 20]; %【受力结点编号  x方向 y方向】 案例为对三角形顶部的1号结点施加竖直向下（y轴方向）的20单位力
FixPos = [22 1 1;28 1 1];%【约束点 x约束 y约束】（有约束为1，无约束为0）

%单元网格内结点编号
NoDsNum =[ 
                                           1 2 3;
                                    2 4 5; 2 5 3; 3 5 6; 
                             4 7 8; 4 8 5; 5 8 9; 5 9 6; 6 9 10;
                    7 11 12; 7 12 8; 8 12 13; 8 13 9;9 13 14; 9 14 10; 10 14 15;
          11 16 17; 11 17 12; 12 17 18; 12 18 13; 13 18 19; 13 19 14; 14 19 20;14 20 15; 15 20 21;
16 22 23; 16 23 17; 17 23 24; 17 24 18; 18 24 25; 18 25 19; 19 25 26; 19 26 20; 20 26 27; 20 27 21;21 27 28;
   ];
%结点坐标数组
NoDsPos = [
                                        0 0;
                                 -1.15 2; 1.15,2; 
                                -2.3 4; 0,4; 2.3,4;
                         -3.45,6; -1.15,6; 1.15,6; 3.45,6;
                         -4.6,8; -2.3,8; 0,8; 2.3,8; 4.6,8;
                -5.75,10; -3.45,10; -1.15,10; 1.15,10; 3.45,10; 5.75,10;
               -6.9,12; -4.6,12; -2.3,12; 0,12; 2.3,12; 4.6,12; 6.9,12;
         -8.05,14; -5.75,14; -3.45,14; -1.15,14; 1.15,14; 3.45,14; 5.75,14; 8.05,14
    ];
%刚度矩阵生成
ASTIF=zeros(2*PointNum,2*PointNum);
for i=1:EleNum
    Y=[1             Poiss                0;
       Poiss           1                  0;
       0               0  ((1-Poiss)/2)*Young/(1-Poiss*Poiss)];
    A = -det([1, NoDsPos(NoDsNum(i,1),1) ,NoDsPos(NoDsNum(i,1),2); 
            1, NoDsPos(NoDsNum(i,2),1), NoDsPos(NoDsNum(i,2),2) 
            1, NoDsPos(NoDsNum(i,3),1), NoDsPos(NoDsNum(i,3),2)])/2;
    for j = 0:2
        b(j+1)=NoDsPos(NoDsNum(i,(rem((j+1),3))+1),2)-NoDsPos(NoDsNum(i,(rem((j+2),3))+1),2);
        c(j+1)=-NoDsPos(NoDsNum(i,(rem((j+1),3))+1),1)+NoDsPos(NoDsNum(i,(rem((j+2),3))+1),1);
    end
    B=[b(1)  0   b(2)   0  b(3)  0;
        0   c(1)  0   c(2)  0   c(3);
        c(1) b(1) c(2) b(2) c(3) b(3)]/(2*A);
    B1(:,:,i)= B;
    S=Y*B;
    ESTIF=B'*S*Thick*A;
    a = NoDsNum(i,:);
    for j = 1:3
        for k=1:3
            ASTIF((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2)= ASTIF((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2)+ESTIF(j*2-1:j*2,k*2-1:k*2);
        end
    end
end
%加入约束
for i = 1:VFixNum
    if FixPos(i,2)==1
        ASTIF(:,(FixPos(i,1)*2-1))=0;
        ASTIF((FixPos(i,1)*2-1),:)=0;
        ASTIF((FixPos(i,1)*2-1),(FixPos(i,1)*2-1))=1;
    end
end
%生成载荷向量
ASLOD(1:2*PointNum)=0;
for i = 1:ForceNum
    ASLOD((ForcePos(i,1)*2-1):ForcePos(i,1)*2)=ForcePos(i,2:3);
end
%求解内力和结点位移量
invASTIF=inv(ASTIF);
ASDISP=ASLOD*invASTIF;
ELEDISP(1:6)=0;
for i=1:EleNum
    i
    for j =1:3
        ELEDISP(j*2-1:j*2) = ASDISP(NoDsNum(i,j)*2-1:NoDsNum(i,j)*2)
    end
    
    STRESS = Y*B1(:,:,i)*ELEDISP'
end