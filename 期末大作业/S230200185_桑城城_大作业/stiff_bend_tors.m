function [keLoc2] = stiff_bend_tors(keLoc,gaussLoc,gaussWts,...
    detJ,invJ,nnode,kt,G,E,Iy,Iz,i)
%%%%% 计算弯曲和扭转刚度 矩阵 
keLoc2 = keLoc;    
%%%%%% gaussLocations = [-1/sqrt(3), 1/sqrt(3)], weight = [1, 1]
for ii = 1:length(gaussLoc)
    xi = gaussLoc(ii,1);
    [~,nDeriv] = shapeFunct_Beam(xi);
    gaussWt = gaussWts(ii,1);

    Xderiv = nDeriv*invJ;
    B = zeros(1,nnode); B(1:nnode) = Xderiv(:);

    %%%%% 局部坐标系中的单元刚度矩阵
    Ke_tx = kt(i,1)*G(i,1)*(B'*B)*detJ*gaussWt;
    Ke_by = E(i,1)*Iy(i,1)*(B'*B)*detJ*gaussWt;
    Ke_bz = E(i,1)*Iz(i,1)*(B'*B)*detJ*gaussWt;

    rxDofL = [4,10];
    ryDofL = [5,11];
    rzDofL = [6,12];

    keLoc2(rxDofL,rxDofL) = keLoc2(rxDofL,rxDofL) + Ke_tx;
    keLoc2(ryDofL,ryDofL) = keLoc2(ryDofL,ryDofL) + Ke_by;
    keLoc2(rzDofL,rzDofL) = keLoc2(rzDofL,rzDofL) + Ke_bz;
end
end

