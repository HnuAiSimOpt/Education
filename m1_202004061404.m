%作者：吕丽娟
%班级：车辆2003班
%时间：2023/3/8
%第一次作业
%用最小二乘法进行线性拟合
%设置数据样本
x =[0 5 10 15 20 25 30 35 40 45 50];
y = [0.2 5.2 9.8 14.5 20.1 25.5 29.7 36 40.1 45 49.9];

[a,b]=Linear_fitting(x,y);

%定义函数求解一次项系数和常数项
function [a,b]=Linear_fitting(x,y)
n=size(x,2);
sum_xy=sum(x.*y);
sum_x=sum(x);
sum_y=sum(y);
sum_xx=sum(x.*x);
%求解系数a,b
a=(n*sum_xy-sum_x*sum_y)/(n*sum_xx-sum_x^2);
b=(sum_xy-sum_xx*sum_y/sum_x)/(sum_x-n*sum_xx/sum_x);
%求解误差
e=(y.*a+b)-y
%%作图 
hold on
xlabel('x');
ylabel('y');
title('OLS');
        y1=a*x+b;
        plot(x,y,'r+',x,y1)
end
