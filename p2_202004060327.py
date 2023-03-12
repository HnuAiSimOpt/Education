# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 12:14:32 2023

@author: admin
"""
#我们利用最小二乘法对数据进行线性拟合，
#生成测试数据值（可以直接由使用者任意输入哦！）
#由公式计算估计的斜率和截距
#画图并评估Linear regression analyse的预测效果
#最后有个运气预测哦！！！快来试试吧！
#1输入测试直线
import numpy as np
import random
import matplotlib.pyplot as plt
k=0
b=0
while True:
    if k==0:
        k=float(input('请输入直线的斜率，例如：4\n'))
    if b==0:
        b=float(input('请输入直线的截距，例如：3\n'))
    if b!=0 and k!=0:
        break
# 2生成测试样本
x=np.linspace(1,100,100)
for i in x:
    y=k*x+b
    y=y+y*0.05*random.random()*(-1)**random.randint(1,2)
print(type(y))
print(y)
#3计算线性拟合的结果
# 3.1计算和式
x_y_multiple=0
x_square=0
x_average=0
y_average=0
for i in range(len(x)):
    x_y_multiple=x_y_multiple+x[i]*y[i]
    x_square=x_square+x[i]**2
    x_average=x_average+x[i]
    y_average=y_average+y[i]

k_estimate=(x_y_multiple-x_average*y_average/len(x))/(x_square-x_average*x_average/len(x))
b_estimate=y_average/len(x)-k_estimate*x_average/len(x)
print('k的预测值为：',k_estimate)
print('b的预测值为：',b_estimate)
#4画图展示预测的结果
for i in x:
    y_estimate=k_estimate*x+b_estimate
# print(y)
# print(y_estimate)
plt.scatter(x,y,color='b',label='original_data')
plt.plot(x,y_estimate,'r',label='linear equation of regression')
plt.legend()
plt.show()

#5计算误差值，由于数据在（0,100），测试集则选择（-25,25）和（75,125）
x_test=list(np.arange(-25,25,1))
for i in range(75,125):
    x_test.append(i)
y_real=[]
y_predict=[]
error=0
y_a=0

for i in range(len(x_test)):
    y_real_value=k*x_test[i]+b
    y_real.append(y_real_value)
    y_predict_value=k_estimate*x_test[i]+b_estimate
    y_predict.append(y_predict_value)
    error=(y_real_value-y_predict_value)**2+error
    y_a=y_a+y_real_value
y_a=y_a/len(x_test)
# print(y_real)
# print(y_predict)
mse=error/len(x_test)
print("均方根误差值：",mse)
#6番外篇-测试你的运气哦！
error_a=0
for i in range(len(x_test)):
    error_a=error_a+ (y_a-(k*x_test[i]+b))**2
R=1-(error/error_a)
print("恭喜你今天运气值为：",R)







