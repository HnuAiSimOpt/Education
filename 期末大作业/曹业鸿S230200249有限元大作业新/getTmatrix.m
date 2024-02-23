function [ T,T24 ] = getTmatrix(x,y,z)
T = zeros(3,3);
x21 = x(2)-x(1);
y21 = y(2)-y(1);
z21 = z(2)-z(1);
x41 = x(4)-x(1);
y41 = y(4)-y(1);
z41 = z(4)-z(1);
L12 = sqrt(x21^2+y21^2+z21^2);
Lz = sqrt((y21*z41-z21*y41)^2+(z21*x41-x21*z41)^2+(x21*y41-y21*x41)^2);
vx = 1/L12*[x21;y21;z21;];
vz = 1/Lz*[y21*z41-z21*y41;z21*x41-x21*z41;x21*y41-y21*x41;];
vy = cross(vz,vx);
T = [vx,vy,vz]';
T1 = [vy,-vx,vz]';
T6 = [T,zeros(3);zeros(3),T1;];
T12 = [T6,zeros(6);zeros(6),T6;];
T24 = [T12,zeros(12);zeros(12),T12;];

end