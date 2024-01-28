function [a,b,c,d,angleA,angleB,angleC,angleD]=angle(x,y)
%计算得到四边形单元的四个内角大小

e1=sqrt((x(2)-x(3))^2+(y(2)-y(3))^2);
e2=sqrt((x(1)-x(4))^2+(y(1)-y(4))^2);

a=sqrt((x(1)-x(2))^2+(y(1)-y(2))^2);
b=sqrt((x(2)-x(4))^2+(y(2)-y(4))^2);
c=sqrt((x(3)-x(1))^2+(y(3)-y(1))^2);
d=sqrt((x(4)-x(3))^2+(y(4)-y(3))^2);

cosA=(a^2+c^2-e1^2)/(2*a*c);
cosB=(a^2+b^2-e2^2)/(2*a*b);
cosC=(b^2+d^2-e1^2)/(2*b*d);
cosD=(c^2+d^2-e2^2)/(2*d*c);

angleA=acos(cosA);
angleB=acos(cosB);
angleC=acos(cosC);
angleD=acos(cosD);
