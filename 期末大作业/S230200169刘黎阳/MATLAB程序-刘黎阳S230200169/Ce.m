function [Ce1,Ce2]=Ce(e_JXY) % C矩阵

% 初始化Ce1和Ce2
Ce1=zeros(9,4);
Ce2=zeros(9,4);

% 高斯积分数据
gauss_points=[0.932469514203152,0.661209386466265,0.238619186083197,-0.932469514203152,-0.661209386466265,-0.238619186083197];
gw=[0.171324492379170,0.360761573048139,0.467913934572691,0.171324492379170,0.360761573048139,0.467913934572691];
kesi=gauss_points;
ita=gauss_points;

% Ce1和Ce的数值积分
for i=1:6
    for j=1:6
        % 压力插值函数
        Fyp= [1/4*(1-kesi(i))*(1-ita(j))
            1/4*(1+kesi(i))*(1-ita(j))
            1/4*(1+kesi(i))*(1+ita(j))
            1/4*(1-kesi(i))*(1+ita(j))]; 

        % 速度插值函数对kesi和ita的导数
        fy_kesi=[1/4*ita(j)*(kesi(i)-1)*(ita(j)-1)+1/4*kesi(i)*ita(j)*(ita(j)-1)
            -ita(j)*kesi(i)*(ita(j)-1)
            1/4*ita(j)*(kesi(i)+1)*(ita(j)-1)+1/4*kesi(i)*ita(j)*(ita(j)-1)
            1/2*(kesi(i)-1)*(1-ita(j)^2)+1/2*kesi(i)*(1-ita(j)^2)
            -2*kesi(i)*(1-ita(j)^2)
            1/2*(kesi(i)+1)*(1-ita(j)^2)+1/2*kesi(i)*(1-ita(j)^2)
            1/4*ita(j)*(kesi(i)-1)*(ita(j)+1)+1/4*kesi(i)*ita(j)*(ita(j)+1)
            -ita(j)*kesi(i)*(ita(j)+1)
            1/4*ita(j)*(kesi(i)+1)*(ita(j)+1)+1/4*kesi(i)*ita(j)*(ita(j)+1)];
        
        fy_ita= [1/4*kesi(i)*(kesi(i)-1)*(ita(j)-1)+1/4*kesi(i)*ita(j)*(kesi(i)-1)
            1/2*(1-kesi(i)^2)*(ita(j)-1)+1/2*ita(j)*(1-kesi(i)^2)
            1/4*kesi(i)*(kesi(i)+1)*(ita(j)-1)+1/4*kesi(i)*ita(j)*(kesi(i)+1)
            -kesi(i)*ita(j)*(kesi(i)-1)
            -2*(1-kesi(i)^2)*ita(j)
            -kesi(i)*ita(j)*(kesi(i)+1)
            1/4*kesi(i)*(kesi(i)-1)*(ita(j)+1)+1/4*kesi(i)*ita(j)*(kesi(i)-1)
            1/2*(1-kesi(i)^2)*(ita(j)+1)+1/2*ita(j)*(1-kesi(i)^2)
            1/4*kesi(i)*(kesi(i)+1)*(ita(j)+1)+1/4*kesi(i)*ita(j)*(kesi(i)+1)];

        % Jacobi相关计算 
        dx_dkesi=fy_kesi'*e_JXY(:,1);
        dx_dita=fy_ita'*e_JXY(:,1);
        dy_dkesi=fy_kesi'*e_JXY(:,2);
        dy_dita=fy_ita'*e_JXY(:,2);

        Jacobi=[dx_dkesi dy_dkesi
                dx_dita  dy_dita];
        AAAA=inv(Jacobi)*[fy_kesi';fy_ita'];
        fy_x=AAAA(1,:)';
        fy_y=AAAA(2,:)';
        det_Jacobi=det(Jacobi);

        % Ce1和Ce2单元方程子块计算
        Ce1=Ce1+gw(i)*gw(j)*fy_x*Fyp'*det_Jacobi;
        Ce2=Ce2+gw(i)*gw(j)*fy_y*Fyp'*det_Jacobi;
    end
end