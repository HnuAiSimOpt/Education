function [B,J] = MatrixB(coord, xi,eta,zeta)
%返回几何矩阵B和雅可比矩阵J
%节点坐标
X = coord(:,1);
Y = coord(:,2);
Z = coord(:,3);

%计算雅克比矩阵
dN1_dxi = -((eta - 1)*(zeta - 1))/8;
dN2_dxi = ((eta - 1)*(zeta - 1))/8;
dN3_dxi = -((eta + 1)*(zeta - 1))/8;
dN4_dxi =((eta + 1)*(zeta - 1))/8;
dN5_dxi =((eta - 1)*(zeta + 1))/8;
dN6_dxi =-((eta - 1)*(zeta + 1))/8;
dN7_dxi =((eta + 1)*(zeta + 1))/8;
dN8_dxi =-((eta + 1)*(zeta + 1))/8;
dN_dxi = [dN1_dxi; dN2_dxi; dN3_dxi; dN4_dxi; dN5_dxi; dN6_dxi; dN7_dxi; dN8_dxi];

 
dN1_deta =-(xi/8 - 1/8)*(zeta - 1);
dN2_deta =(xi/8 + 1/8)*(zeta - 1);
dN3_deta =-(xi/8 + 1/8)*(zeta - 1);
dN4_deta =(xi/8 - 1/8)*(zeta - 1);
dN5_deta =(xi/8 - 1/8)*(zeta + 1);
dN6_deta =-(xi/8 + 1/8)*(zeta + 1);
dN7_deta =(xi/8 + 1/8)*(zeta + 1);
dN8_deta =-(xi/8 - 1/8)*(zeta + 1);
dN_deta = [dN1_deta; dN2_deta; dN3_deta; dN4_deta; dN5_deta; dN6_deta; dN7_deta; dN8_deta];

dN1_dzeta =-(xi/8 - 1/8)*(eta - 1);
dN2_dzeta =(xi/8 + 1/8)*(eta - 1);
dN3_dzeta =-(xi/8 + 1/8)*(eta + 1);
dN4_dzeta =(xi/8 - 1/8)*(eta + 1);
dN5_dzeta =(xi/8 - 1/8)*(eta - 1);
dN6_dzeta =-(xi/8 + 1/8)*(eta - 1);
dN7_dzeta =(xi/8 + 1/8)*(eta + 1);
dN8_dzeta =-(xi/8 - 1/8)*(eta + 1);
dN_dzeta = [dN1_dzeta; dN2_dzeta; dN3_dzeta; dN4_dzeta; dN5_dzeta; dN6_dzeta; dN7_dzeta; dN8_dzeta];


J = zeros(3);
J(1,1) = dN_dxi'*X;
J(1,2) = dN_dxi'*Y;
J(1,3) = dN_dxi'*Z;

J(2,1) = dN_deta'*X;
J(2,2) = dN_deta'*Y;
J(2,3) = dN_deta'*Z;

J(3,1) = dN_dzeta'*X;
J(3,2) = dN_dzeta'*Y;
J(3,3) = dN_dzeta'*Z;

B = zeros(6,24);

for i = 1:8
    dNi=J\[dN_dxi(i); dN_deta(i); dN_dzeta(i)];
    Bi = [dNi(1), 0, 0;...
            0, dNi(2), 0; ...
            0, 0, dNi(3); ...
           dNi(2), dNi(1), 0; ...
           0, dNi(3), dNi(2);...
           dNi(3), 0, dNi(1)];
    B(:,(i-1)*3+1:(i-1)*3+3) = Bi;
end

end


