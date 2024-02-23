function [a,b,c,d,vol]=cal_der_t10(xcoord,ycoord,zcoord)
%--------------------------------------------------------------------------
% Purpose:
%     Calculate A,B,C,D and element volume
%
% Synopsis:
%     [a,b,c,d,vol]=cal_der_t10(x,y,z)
%
% Variable descriptions:
%     A - ai
%     B - bi
%     C - ci
%     D - di
%     vol - element volume
%     xcoord - x-coordinate values
%     ycoord - y-coordinate values 
%     zcoord - z-coordinate values 
%--------------------------------------------------------------------------
x1=xcoord(1); y1=ycoord(1); z1=zcoord(1);
x2=xcoord(2); y2=ycoord(2); z2=zcoord(2);
x3=xcoord(3); y3=ycoord(3); z3=zcoord(3);
x4=xcoord(4); y4=ycoord(4); z4=zcoord(4);

a=zeros(4,1);
b=zeros(4,1);
c=zeros(4,1);
d=zeros(4,1);

a(1)=det([x2 y2 z2;x3 y3 z3;x4 y4 z4]);
a(2)=-det([x1 y1 z1;x3 y3 z3;x4 y4 z4]);
a(3)=det([x1 y1 z1;x2 y2 z2;x4 y4 z4]);
a(4)=-det([x1 y1 z1;x2 y2 z2;x3 y3 z3]);

b(1)=-det([1 y2 z2;1 y3 z3;1 y4 z4]);
b(2)=det([1 y1 z1;1 y3 z3;1 y4 z4]);
b(3)=-det([1 y1 z1;1 y2 z2;1 y4 z4]);
b(4)=det([1 y1 z1;1 y2 z2;1 y3 z3]);

c(1)=det([1 x2 z2;1 x3 z3;1 x4 z4]);
c(2)=-det([1 x1 z1;1 x3 z3;1 x4 z4]);
c(3)=det([1 x1 z1;1 x2 z2;1 x4 z4]);
c(4)=-det([1 x1 z1;1 x2 z2;1 x3 z3]);

d(1)=-det([1 x2 y2;1 x3 y3;1 x4 y4]);
d(2)=det([1 x1 y1;1 x3 y3;1 x4 y4]);
d(3)=-det([1 x1 y1;1 x2 y2;1 x4 y4]);
d(4)=det([1 x1 y1;1 x2 y2;1 x3 y3]);

vol=1/6*det([1 x1 y1 z1;1 x2 y2 z2;1 x3 y3 z3;1 x4 y4 z4]);

return