function k_ele=TriangleElementStiffness(E,miu,t,node_ele)
%TriangleElementStiffness This function returns the element
% stiffness matrix for a Triangle(CST)
% element with modulus of elasticity E,
% Poission's ratio miu, constant thickness t,
% node_ele the node coordinate of element.
% The size of the element stiffness
% matrix is 6 x 6.
%---------node coordinate------
x1=node_ele(1,1);                
y1=node_ele(1,2);
x2=node_ele(2,1);                
y2=node_ele(2,2);
x3=node_ele(3,1);                
y3=node_ele(3,2);
%-------------------------------
A=(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2;   %单元面积
a1=x2*y3-y2*x3;
a2=y1*x3-x1*y3;
a3=x1*y2-y1*x2;
b1=y2-y3;
b2=y3-y1;
b3=y1-y2;
c1=x3-x2;
c2=x1-x3;
c3=x2-x1;
B=1/2/A*[b1 0 b2 0 b3 0;
         0 c1 0 c2 0 c3;
         c1 b1 c2 b2 c3 b3];
D=E/(1-miu^2)*[1 miu 0;
               miu 1 0;
               0 0 (1-miu)/2];     %strss/strain matrix for plane stress
% D=E/(1-miu^2)*[1 miu 0;
%                miu 1 0;
%                0 0 (1-miu)/2];     %strss/strain matrix for plane strain
k_ele=t*A*B'*D*B;    %的单元刚度矩阵
