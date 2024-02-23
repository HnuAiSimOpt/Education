function K = C3D4_K(D, four_nodes_matrix)
%----------------------------------------------%
%function:
%   return B matrix for C3D4 element
%input:
%   D:
%       D matrix for 4-three dimensional analysis
%   four_nodes_matrix:
%       [x1, y1, z1]
%       [x2, y2, z2]
%       [x3, y3, z3]
%       [x4, y4, z4]
%----------------------------------------------%
xyz = zeros(4,4);
xyz(:,2:4) = four_nodes_matrix;
xyz(:,1) = [1;1;1;1];
V = det(xyz)/6;
x1 = xyz(1,2);y1 = xyz(1,3);z1 = xyz(1,4);
x2 = xyz(2,2);y2 = xyz(2,3);z2 = xyz(2,4);
x3 = xyz(3,2);y3 = xyz(3,3);z3 = xyz(3,4);
x4 = xyz(4,2);y4 = xyz(4,3);z4 = xyz(4,4);
if(V <= 0)
    error("The V of C3D4 element is less than 0");
end
mbeta1 = [1 y2 z2;1 y3 z3;1 y4 z4];
mbeta2 = [1 y1 z1;1 y3 z3;1 y4 z4];
mbeta3 = [1 y1 z1;1 y2 z2;1 y4 z4];
mbeta4 = [1 y1 z1;1 y2 z2;1 y3 z3];
mgammal = [1 x2 z2;1 x3 z3;1 x4 z4];
mgamma2 = [1 x1 z1;1 x3 z3;1 x4 z4];
mgamma3 = [1 x1 z1;1 x2 z2;1 x4 z4];
mgamma4 = [1 x1 z1;1 x2 z2;1 x3 z3];
mdelta1 = [1 x2 y2;1 x3 y3;1 x4 y4];
mdelta2 = [1 x1 y1;1 x3 y3;1 x4 y4];
mdelta3 = [1 x1 y1;1 x2 y2;1 x4 y4];
mdelta4 = [1 x1 y1;1 x2 y2;1 x3 y3];
beta1 = -1*det(mbeta1);
beta2 = det(mbeta2);
beta3 = -1*det(mbeta3);
beta4 = det(mbeta4);
gamma1 = det(mgammal);
gamma2 = -1*det(mgamma2);
gamma3 = det(mgamma3);
gamma4 = -1*det(mgamma4);
delta1 = -1*det(mdelta1);
delta2 = det(mdelta2);
delta3 = -1*det(mdelta3);
delta4 = det(mdelta4);
B1 = [beta1 0 0;...
      0 gamma1 0;...
      0 0 delta1;...
      gamma1 beta1 0;...
      0 delta1 gamma1;
      delta1 0 beta1];
B2 = [beta2 0 0;...
      0 gamma2 0;...
      0 0 delta2;...
      gamma2 beta2 0;...
      0 delta2 gamma2;
      delta2 0 beta2];
B3 = [beta3 0 0;...
      0 gamma3 0;...
      0 0 delta3;...
      gamma3 beta3 0;...
      0 delta3 gamma3;
      delta3 0 beta3];
B4 = [beta4 0 0;...
      0 gamma4 0;...
      0 0 delta4;...
      gamma4 beta4 0;...
      0 delta4 gamma4;
      delta4 0 beta4];
B = [B1 B2 B3 B4]/(6*V);
K = V*B'*D*B;
end