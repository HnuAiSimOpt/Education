function [angleA,angleB,angleC,angleD]=cal_angle(x,y)

e1=sqrt((x(2)-x(4))^2+(y(2)-y(4))^2);
e2=sqrt((x(1)-x(3))^2+(y(1)-y(3))^2);

a=sqrt((x(1)-x(2))^2+(y(1)-y(2))^2);
b=sqrt((x(2)-x(3))^2+(y(2)-y(3))^2);
c=sqrt((x(3)-x(4))^2+(y(3)-y(4))^2);
d=sqrt((x(4)-x(1))^2+(y(4)-y(1))^2);

cosA1=(d^2+e2^2-c^2)/(2*d*e2);
cosA2=(a^2+e2^2-b^2)/(2*a*e2);
cosB1=(a^2+e1^2-d^2)/(2*a*e1);
cosB2=(b^2+e1^2-c^2)/(2*b*e1);
cosC1=(b^2+e2^2-a^2)/(2*b*e2);
cosC2=(c^2+e2^2-d^2)/(2*c*e2);
cosD1=(c^2+e1^2-b^2)/(2*c*e1);
cosD2=(d^2+e1^2-a^2)/(2*d*e1);

angleA=acos(cosA1)+acos(cosA2);
angleB=acos(cosB1)+acos(cosB2);
angleC=acos(cosC1)+acos(cosC2);
angleD=acos(cosD1)+acos(cosD2);