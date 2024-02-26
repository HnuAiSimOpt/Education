function [BD,B,D]=BDcalc(s,n,t,Loc,E,NU)
%%  B矩阵求解
%    
% % 
dNsnt=[-(1-n)*(1-t)/8,  -(1-s)*(1-t)/8,  -(1-s)*(1-n)/8; 
    (1-n)*(1-t)/8,   -(1+s)*(1-t)/8,  -(1+s)*(1-n)/8;
    (1+n)*(1-t)/8,   (1+s)*(1-t)/8,  -(1+s)*(1+n)/8;
    -(1+n)*(1-t)/8,   (1-s)*(1-t)/8,  -(1-s)*(1+n)/8;
    -(1-n)*(1+t)/8,  -(1-s)*(1+t)/8,   (1-s)*(1-n)/8;
    (1-n)*(1+t)/8,  -(1+s)*(1+t)/8,   (1+s)*(1-n)/8;
    (1+n)*(1+t)/8,   (1+s)*(1+t)/8,   (1+s)*(1+n)/8;
    -(1+n)*(1+t)/8,  (1-s)*(1+t)/8,    (1-s)*(1+n)/8;];
dNsnt=dNsnt';
J=dNsnt*Loc;
detJ=det(J);

dNxyz=J\dNsnt;
B=zeros(6,24);
for ii=1:8 %计算B矩阵
    Bii=[dNxyz(1,ii) 0 0;0 dNxyz(2,ii) 0;0 0 dNxyz(3,ii);
        dNxyz(2,ii) dNxyz(1,ii) 0;
        0 dNxyz(3,ii) dNxyz(2,ii);
        dNxyz(3,ii) 0 dNxyz(1,ii);];
    B(:,3*(ii-1)+1:3*ii)=Bii;
end
D=[1-NU NU NU 0 0 0;NU 1-NU NU 0 0 0;NU NU 1-NU 0 0 0;0 0 0 0.5-NU 0 0;0 0 0 0 0.5-NU 0;0 0 0 0 0 0.5-NU;];
D=D*(E/((1+NU)*(1-2*NU))); %弹性矩阵

BD=detJ*transpose(B)*D*B;  %BD矩阵
