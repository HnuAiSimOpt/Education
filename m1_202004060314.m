% 车辆2004班 
%鲁昂 
%学号202004060314

%在（0，pi/3)内对y=sinx进行线性拟合

%在y=sinx上取21各点，x取值范围为（0，pi/3)
x = 0:0.05:pi/3
y=sin(x)%得到对应的y值

%求拟合函数y=a0x+a1 的系数a0 a1
R=[x',ones(21,1)]
A=R\y'
a0=A(1,1)
a1=A(2,1)
fprintf('拟合后的函数方程为; y= %8.4f x + %8.4f',a0,a1)
%拟合后的一次函数为y=a0x+a1

%计算残差方和
x1=0:0.1:pi/3 %与x取值范围相同，但步长不同。
y1=a0*x1+a1
ESS=0
for i = 1:1:11
    soloESS=0 
    soloESS=(y1(i)-sin(x1(i)))^2
    ESS=ESS +soloESS
end

%绘制图像
y2=a0*x+a1
figure
plot(x,y,'r',x,y2,'b')

title('拟合后的函数方程')
xlabel('x')
ylabel('y')
