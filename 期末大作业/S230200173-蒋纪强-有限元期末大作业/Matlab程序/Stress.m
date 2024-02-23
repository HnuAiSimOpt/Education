function[stress_g]=Stress(disp,ndof,nd,B,D)
% 计算高斯积分点应力
disp_g=[disp(ndof*(nd(1)-1)+1:ndof*nd(1))
        disp(ndof*(nd(2)-1)+1:ndof*nd(2))
        disp(ndof*(nd(3)-1)+1:ndof*nd(3))
        disp(ndof*(nd(4)-1)+1:ndof*nd(4))
        disp(ndof*(nd(5)-1)+1:ndof*nd(5))
        disp(ndof*(nd(6)-1)+1:ndof*nd(6))
        disp(ndof*(nd(7)-1)+1:ndof*nd(7))
        disp(ndof*(nd(8)-1)+1:ndof*nd(8))];
strain_g=B*disp_g;
stress_g=D*strain_g;