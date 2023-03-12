%% 画出需要拟合的数据点

x=[0.1;0.3;0.4;0.79;0.9;1.1];
y=[1.7805;2.2285;2.3941;3.2226;3.5697;3.8];
plot(x,y,'o')
% xaxis(0,1);
% yaxis(1.6,3.7);
 xlim([0 1.5]);
 ylim([1.6 4]);
hold on

%% 根据最小二乘法原理推导的求解K值和b 进行拟合
% 代数计算
N=length(x);
k=(sum(y.*x)-N*mean(y)*mean(x))/(sum(x.^2)-N*mean(x)^2);
b=mean(y)-k*mean(x);
x1=linspace(0,1.5);
y1=k*x1+b;
plot(x1,y1,'color','r','LineWidth',1.5);


%% 利用自带函数
fun=@(K,x)K(1)*x+K(2);
K0=[1,1];
K=lsqcurvefit(fun,K0,x,y)
k=K(1);
b=K(2);
x_line=linspace(0,1,101);
y_line=k*x_line+b;
plot(x_line,y_line,"Color",'b','LineWidth',1.1);

