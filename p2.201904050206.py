# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 18:34:11 2023
@author: 86186
"""



# 引用库函数

import numpy as np
import random
import matplotlib.pyplot as plt


plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文
plt.rcParams['axes.unicode_minus'] = False    # 用来正常显示负号

# 需要拟合的函数
def f_1(x, A, B,C):
    return A * x**2 + B *x + C

# 需要拟合的数据组
x_group = np.random.randint(0,100,50)
x_group.sort()
y_group = []
y_group = f_1(x_group,1/5,20,-3)
def p_x(x,k,b):    
    return k*x+b
#方程个数
m = len(x_group) 
##最小二乘法求k,b计算过程
sum_x = np.sum(x_group)
sum_y = np.sum(y_group)
sum_xy = np.sum(x_group * y_group)
sum_xx = np.sum(x_group **2 )

b=(sum_y*sum_xx-sum_x*sum_xy)/(m*sum_xx-(sum_x)**2)
k=(m*sum_xy-sum_x*sum_y)/(m*sum_xx-(sum_x)**2)
plt.scatter(x_group, y_group, marker='o',label='真实值')#原始数据散点图
x = x_group
y=p_x(x,k,b)
s=y_group-y
plt.plot(x, y,color='red',label='拟合曲线')#拟合曲线绘图
plt.plot(x,s,color='green',label='误差')
plt.title('y={}*x+{}'.format(k,b)) #拟合曲线表达式
plt.legend() # 显示label
plt.show()
'''由于数据随机生成，故每次运行绘图结果均不同'''
