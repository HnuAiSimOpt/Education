function k = Hexahedral3D8Node_Stiffness(E,NU,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8)
%该函数计算单元的刚度矩阵
%输入弹性模量 E，泊松比 NU
%输入 8 个节点的坐标 x1、y1、z1、x2、y2、z2、x3、y3、z3、x4、y4、z4、x5、y5、z5、x6、y6、z6、x7、y7、z7、x8、y8、z8
%输出单元刚度矩阵 k
% 定义局部坐标系

syms s t n; %创建符号变量，便于之后求导
% 定义形函数 N
N1=(1+s)*(1-t)*(1-n)/8;
N2=(1+s)*(1+t)*(1-n)/8;
N3=(1-s)*(1+t)*(1-n)/8;
N4=(1-s)*(1-t)*(1-n)/8;
N5=(1+s)*(1-t)*(1+n)/8;
N6=(1+s)*(1+t)*(1+n)/8;
N7=(1-s)*(1+t)*(1+n)/8;
N8=(1-s)*(1-t)*(1+n)/8;
% 定义坐标变换
x=N1*x1+N2*x2+N3*x3+N4*x4+N5*x5+N6*x6+N7*x7+N8*x8;
y=N1*y1+N2*y2+N3*y3+N4*y4+N5*y5+N6*y6+N7*y7+N8*y8;
z=N1*z1+N2*z2+N3*z3+N4*z4+N5*z5+N6*z6+N7*z7+N8*z8;
% 定义雅可比矩阵
J=[diff(x,s),diff(y,s),diff(z,s);diff(x,t),diff(y,t),diff(z,t);diff(x,n),diff(y,n),diff(z,n)];
Jdet=det(J);

% 定义B矩阵的系数
a=diff(y,t)*diff(z,n)-diff(z,t)*diff(y,n);
b=diff(y,s)*diff(z,n)-diff(z,s)*diff(y,n);
c=diff(y,s)*diff(z,t)-diff(z,s)*diff(y,t);
d=diff(x,t)*diff(z,n)-diff(z,t)*diff(x,n);
e=diff(x,s)*diff(z,n)-diff(z,s)*diff(x,n);
f=diff(x,s)*diff(z,t)-diff(z,s)*diff(x,t);
g=diff(x,t)*diff(y,n)-diff(y,t)*diff(x,n);
h=diff(x,s)*diff(y,n)-diff(y,s)*diff(x,n);
l=diff(x,s)*diff(y,t)-diff(y,s)*diff(x,t);
% 通过循环计算各个矩阵
Ns=[N1,N2,N3,N4,N5,N6,N7,N8];
Bs=sym(zeros(6,3,8));
for i=1:8
    Bs(:,:,i)=[a*diff(Ns(i),s)-b*diff(Ns(i),t)+c*diff(Ns(i),n),0,0;
        0,-d*diff(Ns(i),s)+e*diff(Ns(i),t)-f*diff(Ns(i),n),0; 0,0,g*diff(Ns(i),s)-h*diff(Ns(i),t)+l*diff(Ns(i),n);
        -d*diff(Ns(i),s)+e*diff(Ns(i),t)-f*diff(Ns(i),n),a*diff(Ns(i),s)-b*diff(Ns(i),t)+c*diff(Ns(i),n),0;
        0,g*diff(Ns(i),s)-h*diff(Ns(i),t)+l*diff(Ns(i),n),-d*diff(Ns(i),s)+e*diff(Ns(i),t)-f*diff(Ns(i),n);
        g*diff(Ns(i),s)-h*diff(Ns(i),t)+l*diff(Ns(i),n),0,a*diff(Ns(i),s)-b*diff(Ns(i),t)+c*diff(Ns(i),n)]/Jdet;
end
% 计算 B 矩阵
B=[Bs(:,:,1),Bs(:,:,2),Bs(:,:,3),Bs(:,:,4),Bs(:,:,5),Bs(:,:,6),Bs(:,:,7),Bs(:,:,8)];

% 计算 D 矩阵
D=(E/((1+NU)*(1-2*NU)))*[1-NU,NU,NU,0,0,0;NU,1-NU,NU,0,0,0;NU,NU,1-NU,0,0,0;
    0,0,0,0.5-NU,0,0;0,0,0,0,0.5-NU,0;0,0,0,0,0,0.5-NU];

% 计算单元刚度矩阵 k
BD=Jdet*transpose(B)*D*B;
z=(int(int(int(BD,n,-1,1),t,-1,1),s,-1,1));

k=double(z);  %转化为双精度类型数据输出
 
