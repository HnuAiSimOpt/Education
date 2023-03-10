# -*- coding: utf-8 -*-
"""
Created on Fri Mar 9 12:54:10 2023

@author: 李世杰
"""
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 10, num=20)    #生成20个在0到20之间的随机数
y = 0.5*x+1.5+np.random.random(20)
plt.scatter(x,y)  #显示随机点的坐标

#设立一系列求解k和b的参数
m = len(x)  
i=0
sum_xy=0
sum_y=0
sum_x=0
sum_xx=0
average_x=0;
average_y=0;
while i<m:
    sum_xy+=x[i]*y[i];
    sum_y+=y[i]
    sum_x+=x[i]
    sum_xx+=x[i]*x[i]
    i+=1
average_x=sum_x/m
average_y=sum_y/m
      
#对k和b的求解过程
k=(m*sum_xy-sum_x*sum_y)/(m*sum_xx-sum_x*sum_x)
b=average_y-average_x*k
print('斜率为:',k)
print('截距为:',b)

#得到线性回归方程并绘制图像
y = x * k + b
plt.plot(x,y,color = 'green')

#误差分析
SE = 0
i=0
for i in range(len(x)):
    SE = SE + (y[i]-k*x[i]-b)**2
    i=i+1

print("误差统计量SE=",SE)