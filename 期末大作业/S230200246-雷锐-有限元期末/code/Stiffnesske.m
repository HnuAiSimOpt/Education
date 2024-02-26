function [KE,B,D]=Stiffnesske(E,NU,Loc)
% %
% 目的：等参单元刚度矩阵
% % 
gsx=[-0.7745966692 0 0.7745966692]; %高斯求积坐标点与系数
gsw=[0.55555555556 0.888888888889 0.55555555556];
D=[1-NU NU NU 0 0 0;NU 1-NU NU 0 0 0;NU NU 1-NU 0 0 0;0 0 0 0.5-NU 0 0;0 0 0 0 0.5-NU 0;0 0 0 0 0 0.5-NU;];
D=D*(E/((1+NU)*(1-2*NU))); %弹性矩阵
KE=zeros(24,24);
for ii=1:3 %三维高斯求积
    sx=gsx(ii);
    sw=gsw(ii);
    for jj=1:3
        nx=gsx(jj);
        nw=gsw(jj);
        for kk=1:3
            tx=gsx(kk);
            tw=gsw(kk);
            [BD,B,D]=BDcalc(sx,nx,tx,Loc,E,NU);
            KE=KE+sw*nw*tw*BD;
        end
    end
end
