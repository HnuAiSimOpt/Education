%% 简介
% 车辆2004 杨志豪 2023/3/8 202004060416 
% 目标：对二次函数y=x*2+x用最小二乘进行一次函数拟合，得到参数以及拟合的直线的函数，并求出其均方误差

%% 初始化
% 建立原始函数与数据
clc,clear;%清空环境变量与命令
yangben = 100*rand([1,1100]);%在0-100的范围内生成1100个随机数样本
x = yangben(:,1:100);%将前100个数据放入列表x，后1000个放入列表x0
x0 = yangben(:,101:1100);
x = sort(x);%升序排列
x0 = sort(x0);
y = x.*x+x;%构造实际函数y=x*2+x，得到100个真实点

%% 最小二乘拟合
% 求解一次拟合函数y=a+bx的a，b
n = length(x);%元素（原始数据）的个数
sum_x = sum(x);%全部x的和
sum_y = sum(y);%全部y的和 
sum_x2 = sum(x.*x);%全部x的平方的和
sum_xy = sum(x.*y);%全部x与对应y的乘积的和
a = (sum_xy*sum_x-sum_y*sum_x2)/(sum_x*sum_x-n*sum_x2);%直接套用PPT13页公式求解
b = (sum_xy*n-sum_y*sum_x)/(sum_x2*n-sum_x*sum_x);%拟合函数为y=a+bx

%% 误差分析
y1 = x0.*x0+x0;%用实际函数得出实际值
y2 = a+b*x0;%用拟合函数得出拟合值
error = sum((y2-y1).*(y2-y1));%计算误差

%% 作图
figure
scatter(x,y,'.','r')%将用来拟合的样本点用红点标出
hold on
plot(x0,y2)%将拟合直线在图上画出
title('最小二乘拟合')%设置标签
xlabel('x')
ylabel('y')
xlim([0 100])%限制范围

%% 打印结果
fprintf('拟合出来直线表达式为:y = %4.2f * x + %4.2f \n',b,a)
fprintf('方差为:%4.2f',error)
