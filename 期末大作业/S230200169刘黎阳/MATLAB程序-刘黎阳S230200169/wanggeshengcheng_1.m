<<<<<<< HEAD
clc
clear
clf

%%  ���򼸺γߴ缰���񻮷ֲ���

H=0.04;   %�����ܸ�
L=0.2;    %�����ܳ�
Nx=30;    %ˮƽ���������������ѡ���ܱ�5��������
Ny=8;    %��ֱ�������������
SS=1.4;   %�����ȣ�����Ϊ��������խ���ֵĸ߶�������ָ߶ȵı�ֵ
theta=10; %����ƽ����ת�Ƕ�

%% �ܵ�Ԫ���ͽ����

E=Nx*Ny;                %�ܵ�Ԫ��
Nz=(2*Nx+1)*(2*Ny+1);   % ���ε�Ԫ�������
Nd=(Nx+1)*(Ny+1);       % ���Ե�Ԫ�������

%%  ��Ԫ���

Dx=L/Nx/2;  %ˮƽ����������
Dy=H/Ny/2;  %��ֱ����������

%% ���ֲ�����

AAA=zeros(Ny*2+1,Nx*2+1);
for i=1:2:2*Nx+1
    AAA(1,i)=(i+1)/2;
end
for i=1:Nx
    AAA(1,2*i)=(Nx+1)*(Ny+1)+i;
end

for i=1:2*Nx+1
    AAA(2,i)=(Nx+1)*(Ny+1)+Nx+i;
end

for j=3:2:2*Ny+1
    for i=1:2:2*Nx+1
        AAA(j,i)=(i+1)/2+(Nx+1)*(j-1)/2;
    end
end

for j=3:2:2*Ny+1
    for i=1:Nx
        AAA(j,2*i)=(Nx+1)*(Ny+1)+(Nx+2*Nx+1)*(j-1)/2+i;
    end
end

for j=4:2:2*Ny
    for i=1:2*Nx+1

        AAA(j,i)=(Nx+1)*(Ny+1)+(j/2)*Nx+(2*Nx+1)*(j/2-1)+i;

    end

end

%%  ��������ϵ��

SL1=Nx/5*2+1;
SL2=Nx/5*2;
SL3=Nx/5*2;
SL4=Nx/5*2;
SL5=Nx/5*2;
SL4=round((2*Nx+1)/5);
SL5=round((2*Nx+1)/5);
SL3=2*Nx+1-SL1-SL2-SL4-SL5;
sl=[ones(1,SL1),(1-(1-SS)/SL2):(SS-1)/SL2:SS,ones(1,SL3)*SS,SS+(1-SS)/SL4:(1-SS)/SL4:1,ones(1,SL5)];

%%  �ı��ζ��ε�ԪJXYV����
for i=1:2*Ny+1
    for j=1:2*Nx+1
        JXYV(AAA(i,j),1)=Dx*(j-1);
        JXYV(AAA(i,j),2)=Dy*(i-1)*sl(j);
    end
end

%%  ����ƽ����ת
for i=1:length(JXYV(:,1))
    R=sqrt((JXYV(i,1)+1)^2+JXYV(i,2)^2);
    theta1=atan(JXYV(i,2)/(JXYV(i,1)+1));
    JXYV(i,1)=R*cos(theta/180*pi+theta1);
    JXYV(i,2)=R*sin(theta/180*pi+theta1);
end

%% �ı��ζ��ε�ԪJMV����
k=0;
for i=1:Ny
    for j=1:Nx
        k=k+1;
        JMV(k,:)=[AAA(2*i-1,2*j-1),AAA(2*i-1,2*j),AAA(2*i-1,2*j+1),...
                  AAA(2*i,2*j-1),AAA(2*i,2*j),AAA(2*i,2*j+1),...
                  AAA(2*i+1,2*j-1),AAA(2*i+1,2*j),AAA(2*i+1,2*j+1),]; 
    end
end

%% �ı������Ե�ԪJMP��JXYP����
JMP=[JMV(:,1),JMV(:,3),JMV(:,9),JMV(:,7)];
JXYP=JXYV([1:Nd],:);

%% BP��������

BP1=AAA(1,:);     %�ױ߶���Ϊ1�ű߽�
BP2=AAA(:,2*Nx+1);%���ڶ���Ϊ2�ű߽�
BP3=AAA(2*Ny+1,:);%�ϱ߶���Ϊ3�ű߽�
BP4=AAA(:,1);     %��ڶ���Ϊ4�ű߽�

%% BE��������

thetax1=pi/2-theta/180*pi;    % 1�ű߽��ⷢ�ַ�����X��н�
thetay1=pi-theta/180*pi;      % 1�ű߽��ⷢ�ַ�����Y��н�
thetax2=theta/180*pi;         % 2�ű߽��ⷢ�ַ�����X��н�
thetay2=pi/2-thetax2;         % 2�ű߽��ⷢ�ַ�����Y��н�
thetax3=pi-pi/2+theta/180*pi; % 3�ű߽��ⷢ�ַ�����X��н�
thetay3=pi-theta/180*pi+pi;   % 3�ű߽��ⷢ�ַ�����Y��н�
thetax4=(180+theta)/180*pi;   % 4�ű߽��ⷢ�ַ�����X��н�
thetay4=pi/2+theta/180*pi;    % 4�ű߽��ⷢ�ַ�����Y��н�
AAA1=ones(Nx,1)*cos(thetax1); % 1�ű߽緽������
AAA2=ones(Nx,1)*cos(thetay1); % 1�ű߽緽������
BE1=[[1:Nx]',ones(size([1:Nx]')),AAA1,AAA2];
BBB1=ones(Ny,1)*cos(thetax2); % 2�ű߽緽������
BBB2=ones(Ny,1)*cos(thetay2); % 2�ű߽緽������
BE2=[[Nx:Nx:Ny*Nx]',2*ones(size([1:Ny]')),BBB1,BBB2];
CCC1=ones(Nx,1)*cos(thetax3); % 3�ű߽緽������
CCC2=ones(Nx,1)*cos(thetay3); % 3�ű߽緽������
BE3=[[Nx*(Ny-1)+1:Nx*Ny]',3*ones(size([1:Nx]')),CCC1,CCC2];
DDD1=ones(Ny,1)*cos(thetax4); % 4�ű߽緽������
DDD2=ones(Ny,1)*cos(thetay4); % 4�ű߽緽������
BE4=[[1:Nx:(Ny-1)*Nx+1]',4*ones(size([1:Ny]')),DDD1,DDD2];

% �����ı���������Ƴ���
rectangle_grid(JMP,JXYV);

%% ���������������洢���

clear Dx Dy  H L Nx Ny i j k
clear theta theta1 R  AAA
clear thetax1 thetax2 thetax3 thetax4
clear thetay1 thetay2 thetay3 thetay4
clear AAA1 AAA2 BBB1 BBB2 CCC1 CCC2 DDD1 DDD2
clear SL1 SL2 SL3 SL4 SL5 SS sl
=======
clc
clear
clf

%%  ���򼸺γߴ缰���񻮷ֲ���

H=0.04;   %�����ܸ�
L=0.2;    %�����ܳ�
Nx=30;    %ˮƽ���������������ѡ���ܱ�5��������
Ny=8;    %��ֱ�������������
SS=1.4;   %�����ȣ�����Ϊ��������խ���ֵĸ߶�������ָ߶ȵı�ֵ
theta=10; %����ƽ����ת�Ƕ�

%% �ܵ�Ԫ���ͽ����

E=Nx*Ny;                %�ܵ�Ԫ��
Nz=(2*Nx+1)*(2*Ny+1);   % ���ε�Ԫ�������
Nd=(Nx+1)*(Ny+1);       % ���Ե�Ԫ�������

%%  ��Ԫ���

Dx=L/Nx/2;  %ˮƽ����������
Dy=H/Ny/2;  %��ֱ����������

%% ���ֲ�����

AAA=zeros(Ny*2+1,Nx*2+1);
for i=1:2:2*Nx+1
    AAA(1,i)=(i+1)/2;
end
for i=1:Nx
    AAA(1,2*i)=(Nx+1)*(Ny+1)+i;
end

for i=1:2*Nx+1
    AAA(2,i)=(Nx+1)*(Ny+1)+Nx+i;
end

for j=3:2:2*Ny+1
    for i=1:2:2*Nx+1
        AAA(j,i)=(i+1)/2+(Nx+1)*(j-1)/2;
    end
end

for j=3:2:2*Ny+1
    for i=1:Nx
        AAA(j,2*i)=(Nx+1)*(Ny+1)+(Nx+2*Nx+1)*(j-1)/2+i;
    end
end

for j=4:2:2*Ny
    for i=1:2*Nx+1

        AAA(j,i)=(Nx+1)*(Ny+1)+(j/2)*Nx+(2*Nx+1)*(j/2-1)+i;

    end

end

%%  ��������ϵ��

SL1=Nx/5*2+1;
SL2=Nx/5*2;
SL3=Nx/5*2;
SL4=Nx/5*2;
SL5=Nx/5*2;
SL4=round((2*Nx+1)/5);
SL5=round((2*Nx+1)/5);
SL3=2*Nx+1-SL1-SL2-SL4-SL5;
sl=[ones(1,SL1),(1-(1-SS)/SL2):(SS-1)/SL2:SS,ones(1,SL3)*SS,SS+(1-SS)/SL4:(1-SS)/SL4:1,ones(1,SL5)];

%%  �ı��ζ��ε�ԪJXYV����
for i=1:2*Ny+1
    for j=1:2*Nx+1
        JXYV(AAA(i,j),1)=Dx*(j-1);
        JXYV(AAA(i,j),2)=Dy*(i-1)*sl(j);
    end
end

%%  ����ƽ����ת
for i=1:length(JXYV(:,1))
    R=sqrt((JXYV(i,1)+1)^2+JXYV(i,2)^2);
    theta1=atan(JXYV(i,2)/(JXYV(i,1)+1));
    JXYV(i,1)=R*cos(theta/180*pi+theta1);
    JXYV(i,2)=R*sin(theta/180*pi+theta1);
end

%% �ı��ζ��ε�ԪJMV����
k=0;
for i=1:Ny
    for j=1:Nx
        k=k+1;
        JMV(k,:)=[AAA(2*i-1,2*j-1),AAA(2*i-1,2*j),AAA(2*i-1,2*j+1),...
                  AAA(2*i,2*j-1),AAA(2*i,2*j),AAA(2*i,2*j+1),...
                  AAA(2*i+1,2*j-1),AAA(2*i+1,2*j),AAA(2*i+1,2*j+1),]; 
    end
end

%% �ı������Ե�ԪJMP��JXYP����
JMP=[JMV(:,1),JMV(:,3),JMV(:,9),JMV(:,7)];
JXYP=JXYV([1:Nd],:);

%% BP��������

BP1=AAA(1,:);     %�ױ߶���Ϊ1�ű߽�
BP2=AAA(:,2*Nx+1);%���ڶ���Ϊ2�ű߽�
BP3=AAA(2*Ny+1,:);%�ϱ߶���Ϊ3�ű߽�
BP4=AAA(:,1);     %��ڶ���Ϊ4�ű߽�

%% BE��������

thetax1=pi/2-theta/180*pi;    % 1�ű߽��ⷢ�ַ�����X��н�
thetay1=pi-theta/180*pi;      % 1�ű߽��ⷢ�ַ�����Y��н�
thetax2=theta/180*pi;         % 2�ű߽��ⷢ�ַ�����X��н�
thetay2=pi/2-thetax2;         % 2�ű߽��ⷢ�ַ�����Y��н�
thetax3=pi-pi/2+theta/180*pi; % 3�ű߽��ⷢ�ַ�����X��н�
thetay3=pi-theta/180*pi+pi;   % 3�ű߽��ⷢ�ַ�����Y��н�
thetax4=(180+theta)/180*pi;   % 4�ű߽��ⷢ�ַ�����X��н�
thetay4=pi/2+theta/180*pi;    % 4�ű߽��ⷢ�ַ�����Y��н�
AAA1=ones(Nx,1)*cos(thetax1); % 1�ű߽緽������
AAA2=ones(Nx,1)*cos(thetay1); % 1�ű߽緽������
BE1=[[1:Nx]',ones(size([1:Nx]')),AAA1,AAA2];
BBB1=ones(Ny,1)*cos(thetax2); % 2�ű߽緽������
BBB2=ones(Ny,1)*cos(thetay2); % 2�ű߽緽������
BE2=[[Nx:Nx:Ny*Nx]',2*ones(size([1:Ny]')),BBB1,BBB2];
CCC1=ones(Nx,1)*cos(thetax3); % 3�ű߽緽������
CCC2=ones(Nx,1)*cos(thetay3); % 3�ű߽緽������
BE3=[[Nx*(Ny-1)+1:Nx*Ny]',3*ones(size([1:Nx]')),CCC1,CCC2];
DDD1=ones(Ny,1)*cos(thetax4); % 4�ű߽緽������
DDD2=ones(Ny,1)*cos(thetay4); % 4�ű߽緽������
BE4=[[1:Nx:(Ny-1)*Nx+1]',4*ones(size([1:Ny]')),DDD1,DDD2];

% �����ı���������Ƴ���
rectangle_grid(JMP,JXYV);

%% ���������������洢���

clear Dx Dy  H L Nx Ny i j k
clear theta theta1 R  AAA
clear thetax1 thetax2 thetax3 thetax4
clear thetay1 thetay2 thetay3 thetay4
clear AAA1 AAA2 BBB1 BBB2 CCC1 CCC2 DDD1 DDD2
clear SL1 SL2 SL3 SL4 SL5 SS sl
>>>>>>> 0c2441c35348cc2cb0d4bfa4f93dc736c7d002fd
save msh_sl