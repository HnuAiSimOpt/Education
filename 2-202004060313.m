%% 简介
%{
作者：张家乐 
学号：202004060313 
日期：2023.3.5
解决的问题：OLS拟合。
本代码文件共分5节：
第一节为初始化，生成了样本数据集；
第二节为OLS拟合，即利用LU分解求解正规方程组，得到拟合函数的两个参数；
第三节计算误差SSE;
第四节绘制拟合函数图像与样本点；
第五节打印结果。
%}
%% 初始化
clc,clear;
x=150*rand([1000,1]);%在0-150的范围内生成服从均匀分布的1000个随机数，作为输入
x=sort(x);%升序排列
y=zeros(1000,1);
for i=x
    y=5*i-3;%由函数y=5*X-3生成相应的真实输出
end
data=[x,y];%数据集建立完毕
len=length(data);%获取数据集长度
train_data=zeros(0.6*len,2);%取前60%个样本点添加-10%--10%的噪声用来拟合
for i=1:0.6*len
    train_data(i,:)=[data(i,1),data(i,2)*(1+0.05*(-1+2*rand(1)))];
end
test_data=data(0.6*len+1:end,:);%后40%个样本点用来计算SSE
%% OLS拟合

%求解正规方程组所用到的系数
m=length(train_data);
sum_x=sum(train_data(:,1));
sum_y=sum(train_data(:,2));
sum_x2=sum(train_data(:,1).^2);
sum_xy=sum(train_data(:,1).*train_data(:,2));

%构造系数矩阵A及向量b
A=[m,sum_x;sum_x,sum_x2];
b=[sum_y;sum_xy];

%对系数矩阵A进行LU分解
U=zeros(2,2);L=zeros(2,2);
for step=1:2%L和U矩阵的求解顺序应为先求U的第一行，再求L的第一列，以此类推
    
    for j=1:2%求U的第step行
        if step==1
            U(step,j)=A(step,j);
        else 
            for k=1:i-1
                he=L(step,k)*U(k,j);
            end
            U(step,j)=A(step,j)-he;
        end
    end
    
    for i=1:2%求L的第step列
        if step==1
            L(i,step)=A(i,step)/U(step,step);
        else 
            for k=1:j-1
                he=L(i,k)*U(k,step);
            end
            L(i,step)=(A(i,step)-he)/U(step,step);
        end
    end
end

%利用Ly=b和Ux=y进行求解x
y1=b(1);y2=b(2)-L(2,1)*y1;
x2=y2/U(2,2);x1=(y1-U(1,2)*x2)/U(1,1);%其中x2为直线斜率，x1为直线截距
%% 计算SSE
  
y_yuce=x2*test_data(:,1)+x1;
SSE=sum((y_yuce-test_data(:,2)).^2);

%% 作图
scatter(train_data(:,1),train_data(:,2),'.','MarkerEdgeColor','#A2142F','MarkerFaceColor','#A2142F');hold on;%拟合数据集绘制出散点图
scatter(test_data(:,1),test_data(:,2),'.','MarkerEdgeColor','#EDB120','MarkerFaceColor','#EDB120');hold on;%测试数据集绘制出散点图
syms x
y=@(x) x2*x+x1;%构建拟合直线的匿名函数
fplot(y,[0,150],'Color','#0072BD','linewidth',2);
grid on;
ylabel('x');
xlabel('y');
legend('用来拟合的样本点','用来测试的样本点','拟合曲线');
title('OLS');
set(get(gca, 'Xlabel'),'Fontname','宋体','FontWeight','bold','Fontsize',20);
set(get(gca, 'Ylabel'),'Fontname','宋体','FontWeight','bold','Fontsize',20);
set(get(gca, 'legend'),'Fontname','宋体','FontWeight','bold','Fontsize',15);
set(get(gca, 'title'),'Fontname','宋体','FontWeight','bold','Fontsize',20);
axis([0 150 0 800]);%为了显示得更清楚，但是坐标轴并不成比例
%% 打印结果
disp(strcat('拟合函数表达式：y=',num2str(x2),'*x',num2str(x1)))
disp('真实函数表达式：y=5*x-3')
disp(strcat('SSE:',num2str(SSE)))
disp('注：每次运行时，样本数据集都会重新随机生成，因此拟合结果有差异')