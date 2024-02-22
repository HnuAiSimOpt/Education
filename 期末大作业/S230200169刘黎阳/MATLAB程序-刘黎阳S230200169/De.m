<<<<<<< HEAD
function [De11,De12,De21,De22]=De(e_JXY,niandu) % D����

%% ��ʼ��De11��De12��De21��De22

De11=zeros(9,9);
De12=zeros(9,9);
De21=zeros(9,9);
De22=zeros(9,9);

%% ��˹��������

gp=[0.932469514203152,0.661209386466265,0.238619186083197,-0.932469514203152,-0.661209386466265,-0.238619186083197];
gw=[0.171324492379170,0.360761573048139,0.467913934572691,0.171324492379170,0.360761573048139,0.467913934572691];
kesi=gp;
ita=gp;

%% De11��De12��De21��De22����ֵ����

for i=1:6
    for j=1:6
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

        %% Jacobi��ؼ��� 

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

        %% De11��De12��De21��De22��Ԫ�����ӿ����

        De11=De11+niandu*gw(i)*gw(j)* [2*fy_x*fy_x'+fy_y*fy_y']*det_Jacobi;
 
        De12=De12+niandu*gw(i)*gw(j)* ...
    fy_x*fy_y'*det_Jacobi;

        De21=De21+niandu*gw(i)*gw(j)* ...
fy_y*fy_x'*det_Jacobi;

        De22=De22+niandu*gw(i)*gw(j)* ...
 [2*fy_y*fy_y'+fy_x*fy_x']*det_Jacobi;

        %%%%%%% De11��De12��De21��De22��Ԫ�����ӿ����

    end

end

=======
function [De11,De12,De21,De22]=De(e_JXY,niandu) % D����

%% ��ʼ��De11��De12��De21��De22

De11=zeros(9,9);
De12=zeros(9,9);
De21=zeros(9,9);
De22=zeros(9,9);

%% ��˹��������

gp=[0.932469514203152,0.661209386466265,0.238619186083197,-0.932469514203152,-0.661209386466265,-0.238619186083197];
gw=[0.171324492379170,0.360761573048139,0.467913934572691,0.171324492379170,0.360761573048139,0.467913934572691];
kesi=gp;
ita=gp;

%% De11��De12��De21��De22����ֵ����

for i=1:6
    for j=1:6
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

        %% Jacobi��ؼ��� 

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

        %% De11��De12��De21��De22��Ԫ�����ӿ����

        De11=De11+niandu*gw(i)*gw(j)* [2*fy_x*fy_x'+fy_y*fy_y']*det_Jacobi;
 
        De12=De12+niandu*gw(i)*gw(j)* ...
    fy_x*fy_y'*det_Jacobi;

        De21=De21+niandu*gw(i)*gw(j)* ...
fy_y*fy_x'*det_Jacobi;

        De22=De22+niandu*gw(i)*gw(j)* ...
 [2*fy_y*fy_y'+fy_x*fy_x']*det_Jacobi;

        %%%%%%% De11��De12��De21��De22��Ԫ�����ӿ����

    end

end

>>>>>>> 0c2441c35348cc2cb0d4bfa4f93dc736c7d002fd
%%%%%%% De11��De12��De21��De22����ֵ����