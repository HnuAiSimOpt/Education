function k_ele=ElementStiffness(E,U,t,node_ele)%6×6单元刚度矩阵

x1=node_ele(1,1);                
y1=node_ele(1,2);
x2=node_ele(2,1);                
y2=node_ele(2,2);
x3=node_ele(3,1);                
y3=node_ele(3,2);

A=(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2;   %三角形单元面积
a1=x2*y3-y2*x3;
a2=y1*x3-x1*y3;
a3=x1*y2-y1*x2;
b1=y2-y3;
b2=y3-y1;
b3=y1-y2;
c1=x3-x2;
c2=x1-x3;
c3=x2-x1;
B=1/2/A.*[b1 0 b2 0 b3 0;0 c1 0 c2 0 c3;c1 b1 c2 b2 c3 b3];%应变矩阵
D=E/(1-U^2).*[1 U 0;U 1 0;0 0 (1-U)/2];   %弹性矩阵
k_ele=t*A*B'*D*B;    %单元刚度矩阵