function [area]=cal_area_3T(xcoord,ycoord,zcoord)
%--------------------------------------------------------------------------
% Purpose:
%     Calculate triangular's area 
%
% Synopsis:
%     [area]=cal_area_3T(xcoord,ycoord,zcoord)
%
% Variable descriptions:
%     area - element area 
%     xcoord - x-coordinate values for the element 
%     ycoord - y-coordinate values for the element 
%     zcoord - z-coordinate values for the element 
%--------------------------------------------------------------------------

x21=xcoord(2)-xcoord(1);
x32=xcoord(3)-xcoord(2);
x13=xcoord(1)-xcoord(3);

y21=ycoord(2)-ycoord(1);
y32=ycoord(3)-ycoord(2);
y13=ycoord(1)-ycoord(3);

z21=zcoord(2)-zcoord(1);
z32=zcoord(3)-zcoord(2);
z13=zcoord(1)-zcoord(3);

a=sqrt(x21^2+y21^2+z21^2);
b=sqrt(x32^2+y32^2+z32^2);
c=sqrt(x13^2+y13^2+z13^2);
s=(a+b+c)/2;
area=sqrt(s*(s-a)*(s-b)*(s-c));

return