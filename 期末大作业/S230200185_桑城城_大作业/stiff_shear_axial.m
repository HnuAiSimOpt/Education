function [keLoc] = stiff_shear_axial(invJ,nnode,A,E,detJ,ks,G,i)
%%%% 计算刚度矩阵
%%%% 用于剪切和轴向变形
%%%%% 初始化局部坐标系中的单元刚度矩阵
keLoc = zeros(12,12); 

%%%%% xi = 0, weight = 2
%%%%% 用于剪切和轴向刚度矩阵
[shape,nDeriv] = shapeFunct_Beam(0);
gaussWt = 2;

Xderiv = nDeriv*invJ;
B = zeros(1,nnode); B(1:nnode) = Xderiv(:);
%%%%% 局部坐标系中的单元刚度矩阵
Ke_ax = A(i,1)*E(i,1)*(B'*B)*detJ*gaussWt;

B_s1 = [B, -shape'];
Ke_s1 = ks(i,1)*G(i,1)*A(i,1)*(B_s1'*B_s1)*detJ*gaussWt;

B_s2 = [B, shape'];
Ke_s2 = ks(i,1)*G(i,1)*A(i,1)*(B_s2'*B_s2)*detJ*gaussWt;

%%%%% 局部刚度矩阵
uDofL = [1,7];
keLoc(uDofL,uDofL) = keLoc(uDofL,uDofL) + Ke_ax;

s1DofL = [2,8, 6,12];
keLoc(s1DofL,s1DofL) = keLoc(s1DofL,s1DofL) + Ke_s1;

s2DofL = [3,9, 5,11];
keLoc(s2DofL,s2DofL) = keLoc(s2DofL,s2DofL) + Ke_s2;

end

