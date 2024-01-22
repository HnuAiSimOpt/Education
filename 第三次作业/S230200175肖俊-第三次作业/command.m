%学号S230200175
%姓名：肖俊
%线性三角形单元
%利用线性三角形单元计算矩形板节点位移和应力
E=210e6    %弹性模量?

NU=0.1      %泊松比

t=0.1           %板厚

k1=LinearTriangleElementStiffness(E,NU,t,0,4,0,2,1,3,1);

k2=LinearTriangleElementStiffness(E,NU,t,0,4,1,3,2,4,1);

k3=LinearTriangleElementStiffness(E,NU,t,2,4,1,3,2,2,1);

k4=LinearTriangleElementStiffness(E,NU,t,1,3,0,2,2,2,1);

k5=LinearTriangleElementStiffness(E,NU,t,2,4,2,2,3,3,1);

k6=LinearTriangleElementStiffness(E,NU,t,2,4,3,3,4,4,1);

k7=LinearTriangleElementStiffness(E,NU,t,3,3,4,2,4,4,1);

k8=LinearTriangleElementStiffness(E,NU,t,3,3,2,2,4,2,1);

k9=LinearTriangleElementStiffness(E,NU,t,2,2,2,0,3,1,1);

k10=LinearTriangleElementStiffness(E,NU,t,2,2,3,1,4,2,1);

k11=LinearTriangleElementStiffness(E,NU,t,3,1,4,0,4,2,1);

k12=LinearTriangleElementStiffness(E,NU,t,3,1,2,0,4,0,1);

k13=LinearTriangleElementStiffness(E,NU,t,4,4,4,2,5,3,1);

k14=LinearTriangleElementStiffness(E,NU,t,4,4,5,3,6,4,1);

k15=LinearTriangleElementStiffness(E,NU,t,5,3,6,2,6,4,1);

k16=LinearTriangleElementStiffness(E,NU,t,4,2,6,2,5,3,1);

K=zeros(28,28); 

K=LinearTriangleAssemble(K,k1,1,4,3);

K=LinearTriangleAssemble(K,k2,1,3,2);

K=LinearTriangleAssemble(K,k3,2,3,5);

K=LinearTriangleAssemble(K,k4,3,4,5);

K=LinearTriangleAssemble(K,k5,2,5,6);

K=LinearTriangleAssemble(K,k6,2,6,8);

K=LinearTriangleAssemble(K,k7,6,7,8);

K=LinearTriangleAssemble(K,k8,6,5,7);

K=LinearTriangleAssemble(K,k9,5,10,9);

K=LinearTriangleAssemble(K,k10,5,9,7);

K=LinearTriangleAssemble(K,k11,9,11,7);

K=LinearTriangleAssemble(K,k12,9,10,11);

K=LinearTriangleAssemble(K,k13,8,7,12);

K=LinearTriangleAssemble(K,k14,8,12,14);

K=LinearTriangleAssemble(K,k15,12,13,14);

K=LinearTriangleAssemble(K,k16,7,13,12)

K(1,1)=K(1,1)*1e15;

K(2,2)=K(2,2)*1e15;

K(19,19)=K(19,19)*1e15;

K(20,20)=K(20,20)*1e15;

K(21,21)=K(21,21)*1e15;

K(22,22)=K(22,22)*1e15;

F=[0;-100;0;-200;0;0;0;0;0;0;0;0;0;0;0;-200;0;0;0;0;0;0;0;0;0;0;0;-100];

U=K\F;

u1=[U(1);U(2);U(7);U(8);U(5);U(6)];

u2=[U(1);U(2);U(5);U(6);U(3);U(4)];

u3=[U(3);U(4);U(5);U(6);U(9);U(10)];

u4=[U(5);U(6);U(7);U(8);U(9);U(10)];

u5=[U(3);U(4);U(9);U(10);U(11);U(12)];

u6=[U(3);U(4);U(11);U(12);U(15);U(16)];

u7=[U(11);U(12);U(13);U(14);U(15);U(16)];

u8=[U(11);U(12);U(9);U(10);U(13);U(14)];

u9=[U(9);U(10);U(19);U(20);U(17);U(18)];

u10=[U(9);U(10);U(17);U(18);U(13);U(14)];

u11=[U(17);U(18);U(21);U(22);U(13);U(14)];

u12=[U(17);U(18);U(19);U(20);U(21);U(22)];

u13=[U(15);U(16);U(13);U(14);U(23);U(24)];

u14=[U(15);U(16);U(23);U(24);U(27);U(28)];

u15=[U(23);U(24);U(25);U(26);U(27);U(28)];

u16=[U(13);U(14);U(25);U(26);U(23);U(24)]

sigma1=LinearTriangleElementStresses(E,NU,0,4,0,2,1,3,1,u1);

sigma2=LinearTriangleElementStresses(E,NU,0,4,1,3,2,4,1,u2);

sigma3=LinearTriangleElementStresses(E,NU,2,4,1,3,2,2,1,u3);

sigma4=LinearTriangleElementStresses(E,NU,1,3,0,2,2,2,1,u4);

sigma5=LinearTriangleElementStresses(E,NU,2,4,2,2,3,3,1,u5);

sigma6=LinearTriangleElementStresses(E,NU,2,4,3,3,4,4,1,u6);

sigma7=LinearTriangleElementStresses(E,NU,3,3,4,2,4,4,1,u7);

sigma8=LinearTriangleElementStresses(E,NU,3,3,2,2,4,2,1,u8);

sigma9=LinearTriangleElementStresses(E,NU,2,2,2,0,3,1,1,u9);

sigma10=LinearTriangleElementStresses(E,NU,2,2,3,1,4,2,1,u10);

sigma11=LinearTriangleElementStresses(E,NU,3,1,4,0,4,2,1,u11);

sigma12=LinearTriangleElementStresses(E,NU,3,1,2,0,4,0,1,u12);

sigma13=LinearTriangleElementStresses(E,NU,4,4,4,2,5,3,1,u13);

sigma14=LinearTriangleElementStresses(E,NU,4,4,5,3,6,4,1,u14);

sigma15=LinearTriangleElementStresses(E,NU,5,3,6,2,6,4,1,u15);

sigma16=LinearTriangleElementStresses(E,NU,4,2,6,2,5,3,1,u16);

s1=LinearTriangleElementPstresses(sigma1);

s2=LinearTriangleElementPstresses(sigma2);

s3=LinearTriangleElementPstresses(sigma3);

s4=LinearTriangleElementPstresses(sigma4);

s5=LinearTriangleElementPstresses(sigma5);

s6=LinearTriangleElementPstresses(sigma6);

s7=LinearTriangleElementPstresses(sigma7);

s8=LinearTriangleElementPstresses(sigma8);

s9=LinearTriangleElementPstresses(sigma9);

s10=LinearTriangleElementPstresses(sigma10);

s11=LinearTriangleElementPstresses(sigma11);

s12=LinearTriangleElementPstresses(sigma12);

s13=LinearTriangleElementPstresses(sigma13);

s14=LinearTriangleElementPstresses(sigma14);

s15=LinearTriangleElementPstresses(sigma15);

s16=LinearTriangleElementPstresses(sigma16);

yipslong1=LinearTriangleElementStresses(E,NU,0,4,0,2,1,3,2,u1);

yipslong2=LinearTriangleElementStresses(E,NU,0,4,1,3,2,4,2,u2);

yipslong3=LinearTriangleElementStresses(E,NU,2,4,1,3,2,2,2,u3);

yipslong4=LinearTriangleElementStresses(E,NU,1,3,0,2,2,2,2,u4);

yipslong5=LinearTriangleElementStresses(E,NU,2,4,2,2,3,3,2,u5);

yipslong6=LinearTriangleElementStresses(E,NU,2,4,3,3,4,4,2,u6);

yipslong7=LinearTriangleElementStresses(E,NU,3,3,4,2,4,4,2,u7);

yipslong8=LinearTriangleElementStresses(E,NU,3,3,2,2,4,2,2,u8);

yipslong9=LinearTriangleElementStresses(E,NU,2,2,2,0,3,1,2,u9);

yipslong10=LinearTriangleElementStresses(E,NU,2,2,3,1,4,2,2,u10);

yipslong11=LinearTriangleElementStresses(E,NU,3,1,4,0,4,2,2,u11);

yipslong12=LinearTriangleElementStresses(E,NU,3,1,2,0,4,0,2,u12);

yipslong13=LinearTriangleElementStresses(E,NU,4,4,4,2,5,3,2,u13);

yipslong14=LinearTriangleElementStresses(E,NU,4,4,5,3,6,4,2,u14);

yipslong15=LinearTriangleElementStresses(E,NU,5,3,6,2,6,4,2,u15);

yipslong16=LinearTriangleElementStresses(E,NU,4,2,6,2,5,3,2,u16);

gElementStress1=[sigma1(1);sigma2(1);sigma3(1);sigma4(1);sigma5(1)

sigma6(1);sigma7(1);sigma8(1);sigma9(1);sigma10(1);sigma11(1)

sigma12(1);sigma13(1);sigma14(1);sigma15(1);sigma16(1)];

gElementStress2=[sigma1(2);sigma2(2);sigma3(2);sigma4(2);sigma5(2)

sigma6(2);sigma7(2);sigma8(2);sigma9(2);sigma10(2);sigma11(2)

sigma12(2);sigma13(2);sigma14(2);sigma15(2);sigma16(2)];

gElementStress3=[sigma1(3);sigma2(3);sigma3(3);sigma4(3)

sigma5(3);sigma6(3);sigma7(3);sigma8(3);sigma9(3);sigma10(3);

sigma11(3);sigma12(3);sigma13(3);sigma14(3);sigma15(3);sigma16(3)];

gElementStress4=[s1(1);s2(1);s3(1);s4(1);s5(1);s6(1);s7(1);s8(1)

s9(1);s10(1);s11(1);s12(1);s13(1);s14(1);s15(1);s16(1)];

gElementStress5=[s1(2);s2(2);s3(2);s4(2);s5(2);s6(2);s7(2);s8(2)

s9(2);s10(2);s11(2);s12(2);s13(2);s14(2);s15(2);s16(2)];

gElementStress=[gElementStress1 gElementStress2 gElementStress3 gElementStress4 gElementStress5];

gElementcoordinate=[0,4,0,2,1,3;0,4,1,3,2,4;2,4,1,3,2,2;1,3,0,2,2,2

2,4,2,2,3,3;2,4,3,3,4,4;3,3,4,2,4,4;3,3,2,2,4,2;2,2,2,0,3,1

2,2,3,1,4,2;3,1,4,0,4,2;3,1,2,0,4,0;4,4,4,2,5,3;4,4,5,3,6,4

5,3,6,2,6,4;4,2,6,2,5,3];