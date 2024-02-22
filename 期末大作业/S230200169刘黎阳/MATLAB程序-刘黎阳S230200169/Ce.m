<<<<<<< HEAD
function [Ce1,Ce2]=Ce(e_JXY) % C����

% ��ʼ��Ce1��Ce2
Ce1=zeros(9,4);
Ce2=zeros(9,4);

% ��˹��������
gauss_points=[0.932469514203152,0.661209386466265,0.238619186083197,-0.932469514203152,-0.661209386466265,-0.238619186083197];
gw=[0.171324492379170,0.360761573048139,0.467913934572691,0.171324492379170,0.360761573048139,0.467913934572691];
kesi=gauss_points;
ita=gauss_points;

% Ce1��Ce����ֵ����
for i=1:6
    for j=1:6
        % ѹ����ֵ����
        Fyp= [1/4*(1-kesi(i))*(1-ita(j))
            1/4*(1+kesi(i))*(1-ita(j))
            1/4*(1+kesi(i))*(1+ita(j))
            1/4*(1-kesi(i))*(1+ita(j))]; 

        % �ٶȲ�ֵ������kesi��ita�ĵ���
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

        % Jacobi��ؼ��� 
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

        % Ce1��Ce2��Ԫ�����ӿ����
        Ce1=Ce1+gw(i)*gw(j)*fy_x*Fyp'*det_Jacobi;
        Ce2=Ce2+gw(i)*gw(j)*fy_y*Fyp'*det_Jacobi;
    end
=======
function [Ce1,Ce2]=Ce(e_JXY) % C����

% ��ʼ��Ce1��Ce2
Ce1=zeros(9,4);
Ce2=zeros(9,4);

% ��˹��������
gauss_points=[0.932469514203152,0.661209386466265,0.238619186083197,-0.932469514203152,-0.661209386466265,-0.238619186083197];
gw=[0.171324492379170,0.360761573048139,0.467913934572691,0.171324492379170,0.360761573048139,0.467913934572691];
kesi=gauss_points;
ita=gauss_points;

% Ce1��Ce����ֵ����
for i=1:6
    for j=1:6
        % ѹ����ֵ����
        Fyp= [1/4*(1-kesi(i))*(1-ita(j))
            1/4*(1+kesi(i))*(1-ita(j))
            1/4*(1+kesi(i))*(1+ita(j))
            1/4*(1-kesi(i))*(1+ita(j))]; 

        % �ٶȲ�ֵ������kesi��ita�ĵ���
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

        % Jacobi��ؼ��� 
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

        % Ce1��Ce2��Ԫ�����ӿ����
        Ce1=Ce1+gw(i)*gw(j)*fy_x*Fyp'*det_Jacobi;
        Ce2=Ce2+gw(i)*gw(j)*fy_y*Fyp'*det_Jacobi;
    end
>>>>>>> 0c2441c35348cc2cb0d4bfa4f93dc736c7d002fd
end