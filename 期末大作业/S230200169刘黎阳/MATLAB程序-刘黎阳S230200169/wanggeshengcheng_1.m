clc
clear
clf

%%  区域几何尺寸及网格划分参数

H=0.04;   %区域总高
L=0.2;    %区域总长
Nx=30;    %水平方向的网格数量，选择能被5整除的数
Ny=8;    %竖直方向的网格数量
SS=1.4;   %收敛比，定义为区域中最窄部分的高度与最宽部分高度的比值
theta=10; %网格平面旋转角度

%% 总单元数和结点数

E=Nx*Ny;                %总单元数
Nz=(2*Nx+1)*(2*Ny+1);   % 二次单元结点总数
Nd=(Nx+1)*(Ny+1);       % 线性单元结点总数

%%  单元间距

Dx=L/Nx/2;  %水平方向网格间距
Dy=H/Ny/2;  %竖直方向网格间距

%% 结点分布拓扑

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

%%  收敛比例系数

SL1=Nx/5*2+1;
SL2=Nx/5*2;
SL3=Nx/5*2;
SL4=Nx/5*2;
SL5=Nx/5*2;
SL4=round((2*Nx+1)/5);
SL5=round((2*Nx+1)/5);
SL3=2*Nx+1-SL1-SL2-SL4-SL5;
sl=[ones(1,SL1),(1-(1-SS)/SL2):(SS-1)/SL2:SS,ones(1,SL3)*SS,SS+(1-SS)/SL4:(1-SS)/SL4:1,ones(1,SL5)];

%%  四边形二次单元JXYV生成
for i=1:2*Ny+1
    for j=1:2*Nx+1
        JXYV(AAA(i,j),1)=Dx*(j-1);
        JXYV(AAA(i,j),2)=Dy*(i-1)*sl(j);
    end
end

%%  网格平面旋转
for i=1:length(JXYV(:,1))
    R=sqrt((JXYV(i,1)+1)^2+JXYV(i,2)^2);
    theta1=atan(JXYV(i,2)/(JXYV(i,1)+1));
    JXYV(i,1)=R*cos(theta/180*pi+theta1);
    JXYV(i,2)=R*sin(theta/180*pi+theta1);
end

%% 四边形二次单元JMV生成
k=0;
for i=1:Ny
    for j=1:Nx
        k=k+1;
        JMV(k,:)=[AAA(2*i-1,2*j-1),AAA(2*i-1,2*j),AAA(2*i-1,2*j+1),...
                  AAA(2*i,2*j-1),AAA(2*i,2*j),AAA(2*i,2*j+1),...
                  AAA(2*i+1,2*j-1),AAA(2*i+1,2*j),AAA(2*i+1,2*j+1),]; 
    end
end

%% 四边形线性单元JMP和JXYP生成
JMP=[JMV(:,1),JMV(:,3),JMV(:,9),JMV(:,7)];
JXYP=JXYV([1:Nd],:);

%% BP数据生成

BP1=AAA(1,:);     %底边定义为1号边界
BP2=AAA(:,2*Nx+1);%出口定义为2号边界
BP3=AAA(2*Ny+1,:);%上边定义为3号边界
BP4=AAA(:,1);     %入口定义为4号边界

%% BE数据生成

thetax1=pi/2-theta/180*pi;    % 1号边界外发现方向与X轴夹角
thetay1=pi-theta/180*pi;      % 1号边界外发现方向与Y轴夹角
thetax2=theta/180*pi;         % 2号边界外发现方向与X轴夹角
thetay2=pi/2-thetax2;         % 2号边界外发现方向与Y轴夹角
thetax3=pi-pi/2+theta/180*pi; % 3号边界外发现方向与X轴夹角
thetay3=pi-theta/180*pi+pi;   % 3号边界外发现方向与Y轴夹角
thetax4=(180+theta)/180*pi;   % 4号边界外发现方向与X轴夹角
thetay4=pi/2+theta/180*pi;    % 4号边界外发现方向与Y轴夹角
AAA1=ones(Nx,1)*cos(thetax1); % 1号边界方向余弦
AAA2=ones(Nx,1)*cos(thetay1); % 1号边界方向余弦
BE1=[[1:Nx]',ones(size([1:Nx]')),AAA1,AAA2];
BBB1=ones(Ny,1)*cos(thetax2); % 2号边界方向余弦
BBB2=ones(Ny,1)*cos(thetay2); % 2号边界方向余弦
BE2=[[Nx:Nx:Ny*Nx]',2*ones(size([1:Ny]')),BBB1,BBB2];
CCC1=ones(Nx,1)*cos(thetax3); % 3号边界方向余弦
CCC2=ones(Nx,1)*cos(thetay3); % 3号边界方向余弦
BE3=[[Nx*(Ny-1)+1:Nx*Ny]',3*ones(size([1:Nx]')),CCC1,CCC2];
DDD1=ones(Ny,1)*cos(thetax4); % 4号边界方向余弦
DDD2=ones(Ny,1)*cos(thetay4); % 4号边界方向余弦
BE4=[[1:Nx:(Ny-1)*Nx+1]',4*ones(size([1:Ny]')),DDD1,DDD2];

% 调用四边形网格绘制程序
rectangle_grid(JMP,JXYV);

%% 清除多余变量，并存储结果

clear Dx Dy  H L Nx Ny i j k
clear theta theta1 R  AAA
clear thetax1 thetax2 thetax3 thetax4
clear thetay1 thetay2 thetay3 thetay4
clear AAA1 AAA2 BBB1 BBB2 CCC1 CCC2 DDD1 DDD2
clear SL1 SL2 SL3 SL4 SL5 SS sl
save msh_sl