function y = LinearTriangleElementPStresses(sigma, p)
%LinearTriangleElementPStresses    函数给出单元主应力
%                                  与主应力方向
%                                  给定应力张量

R = (sigma(1)+sigma(2))/2;
Q = ((sigma(1)-sigma(2))/2)^2+sigma(3)*sigma(3);
M = 2*sigma(3)/(sigma(1)-sigma(2));
s1 = R+sqrt(Q);
s2 = R-sqrt(Q);
theta = (atan(M)/2)*180/pi;
if p ==1
    y = [s1;s2;theta];
elseif p==2
    y = [1/sqrt(2)*sqrt((s1-s2)^2+s1^2+s2^2), 0, 0];
end