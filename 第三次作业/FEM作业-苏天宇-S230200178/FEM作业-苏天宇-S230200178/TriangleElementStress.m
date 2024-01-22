function str=TriangleElementStress(E,u,node_ele,u1,p)
x1=node_ele(1,1);                
y1=node_ele(1,2);
x2=node_ele(2,1);                
y2=node_ele(2,2);
x3=node_ele(3,1);                
y3=node_ele(3,2);
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
if p==1
    D=E/(1-u^2)*[1 u 0;
                   u 1 0;
                   0 0 (1-u)/2];           
elseif p==2
    D=E/(1+u)/(1-2*u)*[1-u u 0;
                       u 1-u 0;
                       0 0 (1-2*u)/2];   
end
str=D*B*u1;   %单元应力
end