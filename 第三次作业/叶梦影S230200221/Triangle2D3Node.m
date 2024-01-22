%姓名：叶梦影
%学号：S230200221
%一矩形薄平板长2m，宽1m，在右端部受集中力F=100KN，弹性模量E=1e7Pa，泊松比1/3，板的厚度t=0.1m。
%求解该结构的节点位移、支反力以及单元应力。(将该矩形单元分成两个3节点三角形单元)
%----------------------------------------------------------
E=1e7;
NU=1/3;
t=0.1;
ID=1;
k1=Triangle2D3Node_Stiffness(E,NU,t,2,0,0,1,0,0,ID);
k2=Triangle2D3Node_Stiffness(E,NU,t,0,1,2,0,2,1,ID);
KK = zeros(8,8);
KK=Triangle2D3Node_Assembly(KK,k1,2,3,4);
KK=Triangle2D3Node_Assembly(KK,k2,3,2,1);
k=KK(1:4,1:4);
%针对1，2节点的位移进行求解，3，4节点两个方向上的位移为0
p=[0;-5000;0;-5000];
u=k\p;
fprintf('节点位移为u=\n');
disp(u);
%求支反力-------------------------
U=[u;0;0;0;0];
P=KK*U;
fprintf('支反力为：P=\n');
disp(P);%对应关系P=[0; -F/2; 0; -F/2; Rx3; Ry3; Rx4; Ry4]
%各单元应力的计算---------------------------------------
u1=[U(3);U(4);U(5);U(6);U(7);U(8)];
stress1=Triangle2D3Node_Stress(E,NU,2,0,0,1,0,0,u1,ID);
%输出单元1的应力 stress(3X1)，由于它为常应力单元，则单元的应力分量为 Sx,Sy,Sxy 
fprintf('单元1的应力为：stress1=\n');
disp(stress1);
u2=[U(5);U(6);U(3);U(4);U(1);U(2)];
stress2=Triangle2D3Node_Stress(E,NU,0,1,2,0,2,1,u2,ID);
fprintf('单元2的应力为：stress1=\n');
disp(stress2);
%输出单元2的应力 stress(3X1)，由于它为常应力单元，则单元的应力分量为 Sx,Sy,Sxy 

function k=Triangle2D3Node_Stiffness(E,NU,t,xi,yi,xj,yj,xm,ym,ID) 
%该函数计算单元的刚度矩阵 
%输入弹性模量 E，泊松比 NU，厚度 t 
%输入三个节点 i、j、m 的坐标 xi,yi,xj,yj,xm,ym  
A = (xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2; 
betai = yj-ym; 
betaj = ym-yi; 
betam = yi-yj; 
gammai = xm-xj; 
gammaj = xi-xm; 
gammam = xj-xi; 
B = [betai 0 betaj 0 betam 0 ; 
0 gammai 0 gammaj 0 gammam ; 
gammai betai gammaj betaj gammam betam]/(2*A); 
if ID == 1 
D = (E/(1-NU*NU))*[1 NU 0 ; NU 1 0 ; 0 0 (1-NU)/2]; 
elseif ID == 2 
D = (E/(1+NU)/(1-2*NU))*[1-NU NU 0 ; NU 1-NU 0 ; 0 0 (1-2*NU)/2]; 
end 
k= t*A*B'*D*B; 
end
%--------------------------------------------------------------- 
function z = Triangle2D3Node_Assembly(KK,k,i,j,m) 
%该函数进行单元刚度矩阵的组装 
%输入单元刚度矩阵 k 
%输入单元的节点编号 I、j、m 
DOF(1)=2*i-1; 
DOF(2)=2*i; 
DOF(3)=2*j-1; 
DOF(4)=2*j; 
DOF(5)=2*m-1; 
DOF(6)=2*m; 
for n1=1:6 
for n2=1:6 
KK(DOF(n1),DOF(n2))= KK(DOF(n1),DOF(n2))+k(n1,n2); 
end 
end 
z=KK; 
end
%--------------------------------------------------------------- 
function stress=Triangle2D3Node_Stress(E,NU,xi,yi,xj,yj,xm,ym,u,ID) 
%该函数计算单元的应力 
%输入弹性模量 E，泊松比 NU，厚度 t 
%输入三个节点 i、j、m 的坐标 xi,yi,xj,yj,xm,ym 
%输入平面问题性质指示参数 ID(1 为平面应力，2 为平面应变)，单元的位移列阵 u(6X1) 
A = (xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2; 
betai = yj-ym; 
betaj = ym-yi; 
betam = yi-yj; 
gammai = xm-xj; 
gammaj = xi-xm; 
gammam = xj-xi; 
B = [betai 0 betaj 0 betam 0 ; 
0 gammai 0 gammaj 0 gammam ; 
gammai betai gammaj betaj gammam betam]/(2*A); 
if ID == 1 
D = (E/(1-NU*NU))*[1 NU 0 ; NU 1 0 ; 0 0 (1-NU)/2]; 
elseif ID == 2 
D = (E/(1+NU)/(1-2*NU))*[1-NU NU 0 ; NU 1-NU 0 ; 0 0 (1-2*NU)/2]; 
end 
stress = D*B*u; 
end

