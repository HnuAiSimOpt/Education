function stress=T_Stress(E,NU,xe,ye,ele_u,ID)
%该函数计算单元的应力

if ID == 1 
   D = (E/(1-NU*NU))*[1 NU 0 ; NU 1 0 ; 0 0 (1-NU)/2];
elseif ID == 2
   D = (E/(1+NU)/(1-2*NU))*[1-NU NU 0 ; NU 1-NU 0 ; 0 0 (1-2*NU)/2];
end
xy = zeros(3,2);
for i = 1:3
    xy(i,1)=xe(i);
    xy(i,2)=ye(i);
end
[J B Ae]=T_B(xy);
stress = D*B*ele_u;
