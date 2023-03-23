%% 简介
%{
作者：曾婷曦
学号：202004061326 
日期：2023.3.8
解决的问题：OLS拟合。
本代码文件共分4步：
第一步生成样本数据；
第二步OLS拟合，即利用法方程的矩阵形式得到拟合函数的两个参数；
第三步计算误差SSE;
第四步绘制拟合函数图像与样本点；
%}
%% 初始化
clc,clear;
%数据样本
x=[1,1.5,2,2.5,3.5,4,4.5]; 
y=[0.9,1.7,2.2,2.6,6,10,12];      


%% OLS拟合

n = length(x);%变量个数
sumx = 0;
sumx2 = 0;
sumy = 0;
sumxy = 0;

for i = 1:n
    sumx = sumx + x(i);% 求Σ(xi)
    sumx2 = sumx2 + x(i)*x(i);% 求Σ(xi.^2)
    sumy = sumy + y(i);%求Σ(yi)
    sumxy = sumxy + x(i)*y(i);%求Σ(xi*yi)

end
%构造系数矩阵A及向量b
A=[n,sumx;sumx,sumx2];
b=[sumy;sumxy];

% 求解最小二乘问题Ax=b，得到系数C=[k,b]
C=A\b;
%% 计算SSE
% 拟合直线上的点
xfit =[1,1.5,2,2.5,3.5,4,4.5];
yfit = C(2)*xfit + C(1);
SSE=sum((yfit-y).^2);
%输出拟合公式
disp(['拟合直线：y=' num2str(C(2)) 'x+' num2str(C(1))])
disp(['残差平方和：' SSE])

%% 绘制散点图和拟合直线
figure
plot(x,y,'.')%拟合数据集绘制出散点图
hold on
plot(xfit,yfit,'r-')
xlabel('x')
ylabel('y')
legend('原始数据','拟合直线')
title('OLS');
set(get(gca, 'Xlabel'),'Fontname','宋体','FontWeight','bold','Fontsize',20);
set(get(gca, 'Ylabel'),'Fontname','宋体','FontWeight','bold','Fontsize',20);
set(get(gca, 'legend'),'Fontname','宋体','FontWeight','bold','Fontsize',15);
set(get(gca, 'title'),'Fontname','宋体','FontWeight','bold','Fontsize',20);
