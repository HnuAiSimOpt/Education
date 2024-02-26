function Ke = stiffnessMatrix(coord,E,nu)
%输入单元8个节点坐标coord，弹性模量E和泊松比nu

%采用3*3高斯积分计算单元刚度矩阵
t = [-0.7745966692, 0, 0.7745966692];
w = [0.5555555556,0.8888888889, 0.5555555556];

Ke = zeros(24);
D = MatrixD(E,nu);

for i = 1:length(t)
    for j = 1:length(t)
        for k = 1:length(t)
            xi = t(i);
            eta = t(j);
            zeta = t(k);
            [B,J] = MatrixB(coord,xi,eta,zeta);
            Ke = Ke + w(i)*w(j)*w(k)*transpose(B)*D*B*det(J);
            
        end
    end
end

