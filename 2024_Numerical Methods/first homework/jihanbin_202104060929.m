clear
clc
close all
% 使用2000年苏必利尔湖湖全年水位（按月采样）进行拉格朗日插值
water_level_data = [173.82	173.75	173.83	173.95	174.07	174.18	174.26	174.22	174.09	174	173.88	173.83];
%创建时间坐标（1-12月）
t = 1:12;
t_continue = 1:0.001:12;
syms x
%创建多项式函数l,用x/x使L初始化为可以接收符号变量的行向量
l = x/x*ones(1,12);
%舒适化插值函数L
L=0;
for i =1:12
    for j =1:12
        if i~=j
            l(i)=l(i)*(x-t(j))/(t(i)-t(j));
        end
    end
    L=L+l(i)*water_level_data(i);
end
L = collect(L, x);
%使用matlab自带的多项式函数进行插值作为真值，比较结果
L_t_function = interp1(t,water_level_data,t_continue,'spline');

%可视化
disp("插值函数为：")
disp(L)
L_t_self=subs(L,x,t_continue);
figure(1)
subplot(1,2,1)
plot(t,water_level_data,'or',"DisplayName","观测水位")
hold on 
plot(t_continue,L_t_self,'-b',"DisplayName","插值函数")
xlabel("t /month")
ylabel('waterlevel /mm')
title("手写拉格朗日插值")
legend 

subplot(1,2,2)
plot(t,water_level_data,'or',"DisplayName","观测水位")
hold on 
plot(t_continue,L_t_function,'-b',"DisplayName","插值函数")
xlabel("t /month")
ylabel('waterlevel /mm')
title("matlab多项式插值")
legend 