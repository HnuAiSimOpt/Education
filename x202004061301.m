%许少辉 
% 车辆2003班
% 202004061301 
%2023.03.06
%采用最小二乘算法，通过梯度下降来求损失函数的最小值
%另上传了C语言的代码文件，但用C写的拟合误差太大且找不到问题所在


x=[0 10 20 30 40 50 60 70 80];  %数据样本
y=[10.6 21 29.6 40.5 50 60.7 69.3 80 90.2];
theta1=0;   %一次项系数
theta0=0;   %常系数
v=0.00001; %设置步长
[row,col]=size(x);   %一次项系数和常系数初始化
theta1 = mean(y)/mean(x)
theta0 = y(1)-theta1*x(1)
while 1     %梯度下降，进行循环直到找到最优解退出循环
    temp1=0;%
    temp0=0;
    for i=1:col%对损失函数求导
        temp1 = temp1 - ((y(i) - (theta1*x(i) + theta0)) * x(i));
        temp0=temp0-((y(i)-(theta1*x(i)+theta0))*1);
        end
    old_theta1=theta1;%前一个常系数和一次项系数存储以后续比较
    old_theta0=theta0;
    theta1 = theta1- v*temp1%更新每个样本的常系数和一次项系数
    theta0 =theta0 - v*temp0
    temp1=0;
    temp0=0;
    e =((old_theta1-theta1)^2+(old_theta0-theta0)^2)%误差判别
    if e<0.000003
        f=theta1*x+theta0;%绘图
        plot(x,y,'r+',x,f)
        break;
     end
end