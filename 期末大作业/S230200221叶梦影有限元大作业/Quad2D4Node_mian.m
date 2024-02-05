%***********************************************************************
%       4节点四边形单元主程序
%       姓名：叶梦影
%       学号：S230200221
%***********************************************************************
%自定义外力载荷q0=100KN/m、弹性模量E=7e10、泊松比Nu=0.33、厚度t=0.1m
clc,clear
q0 = 1000000;
E = 7e10;
NU = 0.33;
t = 0.1;

%  输入所划分单元的各个节点的坐标，共25个节点，默认每个单元长度为1m
x1 = 0; y1 = 0;
x2 = 0; y2 = 1;
x3 = 0; y3 = 2;
x4 = 0; y4 = 3;
x5 = 0; y5 = 4;
x6 = 1; y6 = 4;
x7 = 1; y7 = 3;
x8 = 1; y8 = 2;
x9 = 1; y9 = 1;
x10 = 1; y10 = 0;
x11 = 2; y11 = 0;
x12 = 2; y12 = 1;
x13 = 2; y13 = 2;
x14 = 2; y14 = 3;
x15 = 2; y15 = 4;
x16 = 3; y16 = 4;
x17 = 3; y17 = 3;
x18 = 3; y18 = 2;
x19 = 3; y19 = 1;
x20 = 3; y20 = 0;
x21 = 4; y21 = 0;
x22 = 4; y22 = 1;
x23 = 4; y23 = 2;
x24 = 4; y24 = 3;
x25 = 4; y25 = 4;


%计算各个单元的刚度矩阵，共有16个单元
k1 = Quad2D4Node_Stiffness(E, NU, t, x1, y1, x10, y10, x9, y9, x2, y2);
k2 = Quad2D4Node_Stiffness(E, NU, t, x2, y2, x9, y9, x8, y8, x3, y3);
k3 = Quad2D4Node_Stiffness(E, NU, t, x3, y3, x8, y8, x7, y7, x4, y4);
k4 = Quad2D4Node_Stiffness(E, NU, t, x4, y4, x7, y7, x6, y6, x5, y5);
k5 = Quad2D4Node_Stiffness(E, NU, t, x7, y7, x14, y14, x15, y15, x6, y6);
k6 = Quad2D4Node_Stiffness(E, NU, t, x8, y8, x13, y13, x14, y14, x7, y7);
k7 = Quad2D4Node_Stiffness(E, NU, t, x9, y9, x12, y12, x13, y13, x8, y8);
k8 = Quad2D4Node_Stiffness(E, NU, t, x10, y10, x11, y11, x12, y12, x9, y9);
k9 = Quad2D4Node_Stiffness(E, NU, t, x11, y11, x20, y20, x19, y19, x12, y12);
k10 = Quad2D4Node_Stiffness(E, NU, t, x12, y12, x19, y19, x18, y18, x13, y13);
k11 = Quad2D4Node_Stiffness(E, NU, t, x13, y13, x18, y18, x17, y17, x14, y14);
k12 = Quad2D4Node_Stiffness(E, NU, t, x14, y14, x17, y17, x16, y16, x15, y15);
k13 = Quad2D4Node_Stiffness(E, NU, t, x17, y17, x24, y24, x25, y25, x16, y16);
k14 = Quad2D4Node_Stiffness(E, NU, t, x18, y18, x23, y23, x24, y24, x17, y17);
k15 = Quad2D4Node_Stiffness(E, NU, t, x19, y19, x22, y22, x23, y23, x18, y18);
k16 = Quad2D4Node_Stiffness(E, NU, t, x20, y20, x21, y21, x22, y22, x19, y19);

%  组装整体的刚度矩阵
KK = zeros(50,50);

KK = Quad2D4Node_Assembly(KK, k1, 1, 10, 9, 2);
KK = Quad2D4Node_Assembly(KK, k2, 2, 9, 8, 3);
KK = Quad2D4Node_Assembly(KK, k3, 3, 8, 7, 4);
KK = Quad2D4Node_Assembly(KK, k4, 4, 7, 6, 5);
KK = Quad2D4Node_Assembly(KK, k5, 7, 14, 15, 6);
KK = Quad2D4Node_Assembly(KK, k6, 8, 13, 14, 7);
KK = Quad2D4Node_Assembly(KK, k7, 9, 12,13,8);
KK = Quad2D4Node_Assembly(KK, k8, 10, 11,12,9);
KK = Quad2D4Node_Assembly(KK, k9, 11,20,19,12);
KK = Quad2D4Node_Assembly(KK, k10, 12,19,18,13);
KK = Quad2D4Node_Assembly(KK, k11, 13,18,17,14);
KK = Quad2D4Node_Assembly(KK, k12, 14,17,16,15);
KK = Quad2D4Node_Assembly(KK, k13, 17,24,25,16);
KK = Quad2D4Node_Assembly(KK, k14, 18,23,24,17);
KK = Quad2D4Node_Assembly(KK, k15, 19,22,23,18);
KK = Quad2D4Node_Assembly(KK, k16, 20,21,22,19);


% 删除总刚度矩阵中固定约束处，所对应节点位移为0的行和列数据
k=KK(11:50,11:50);

%输入外界所施加力，并且维度匹配，元素个数共有72-13+1=60
p=[ 0;0;0;0;0;0;0;0;0;0;
    0;0;0;0;0;0;0;0;0;0;
    0;0;0;0;0;0;0;0;0;0;
    50000;0;100000;0;100000;0;100000;0;50000;0;];
u = k\p         %高斯消去法求解并输出各节点位移，其中前5个节点位移为0，因此在结果中并不显示

%维度匹配
U = [0;0;0;0;0;0;0;0;0;0;u];

P = KK*U        %输出的P为各节点的载荷，共25个节点，一共输出2*36=50个数据
