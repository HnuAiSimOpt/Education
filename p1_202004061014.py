# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 18:34:11 2023

@author: 86186
"""



# 引用库函数

import numpy as np
import random
import matplotlib.pyplot as plt
from scipy import optimize as op

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文
plt.rcParams['axes.unicode_minus'] = False    # 用来正常显示负号

# 需要拟合的函数
def f_1(x, A, B,C):
    return A * x**2 + B *x + C

# 需要拟合的数据组
x_group = np.random.randint(0,100,50)
x_group.sort()
y_group = np.random.randint(101,500,50)
y_group.sort()


# 得到返回的A，B值
A, B,C = op.curve_fit(f_1, x_group, y_group)[0]

# 数据点与原先的进行画图比较
plt.scatter(x_group, y_group, marker='o',label='真实值')
x = x_group
y = A * x**2 + B *x + C
e=y-y_group#
plt.plot(x, y,color='red',label='拟合曲线')
plt.plot(e,color="green",label="误差")
plt.plot()
plt.title('y={}x^2+{}x+{}'.format(A,B,C))
plt.legend() # 显示label

plt.show()
'''由于数据随机生成，故每次运行绘图结果均不同'''