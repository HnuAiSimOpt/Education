function y =LinearTriangleElementPstresses(sigma)

R=(sigma(1)+sigma(2))/2;

Q=((sigma(1)-sigma(2))/2)^2+sigma(3)*sigma(3);

M=2*sigma(3)/(sigma(1)-sigma(2));

s1=R+sqrt(Q);

s2=R-sqrt(Q);

theta=(atan(M)/2)*180/pi;

y=[s1;s2;theta];
