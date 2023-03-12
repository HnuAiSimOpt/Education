# -*- coding: utf-8 -*-
"""
王洪苑
车辆2003
202014010709
2023.03.07

使用了OLS最小二乘法进行对符合高斯分布数据样本的线性拟合和误差分析
直接按照最小二乘法解方程进行求解
调用了一些库进行高斯分布数据的生成和你和结果的图像展示

@author: Hom Y
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve

# 生成随机25个横坐标 区间在(0,10) 并且转换为矩阵
X = np.linspace(0,10,num = 25).reshape(-1,1)
    
# 随机生成斜率与截距
W = np.random.randint(1,5,size = 1)
B = np.random.randint(1,10,size = 1)

# 根据一元一次方程计算目标值y，加上“噪声”
# 构建初始模型并借此生成符合高斯分布的样本数据，
Y = X * W + B + np.random.randn(25,1)

# 显示随机坐标点
plt.scatter(X,Y)

# 构建求解OLS线性回归系数方程组系数
a1 = 0
for i1 in range(len(X)):
    a1 = a1 + 2*(X[[i1]])**2
b1 = 0
for i2 in range(len(X)):
    b1 = b1 + 2*(X[[i2]])
c1 = 0
for i3 in range(len(X)):
    c1 = c1 + 2*Y[[i3]]*X[[i3]]
a2 = 0
for k1 in range(len(X)):
    a2 = a2 + 2*(X[[k1]])
b2 = 0
for k2 in range(len(X)):
    b2 = b2 + np.array([[2]])
c2 = 0
for k3 in range(len(X)):
    c2 = c2 + 2*Y[[k3]]

#求解方程得到斜率与截距
a=np.mat([[float(a1),float(b1)],[float(a2),float(b2)]])#系数矩阵
b=np.mat([float(c1),float(c2)]).T    #常数项列矩阵
x=solve(a,b)        #方程组的解
print("线性回归方程斜率为：",x[0])
print("           截距为：",x[1])

#得到线性回归方程
y = X * x[0] + x[1]
plt.plot(X,y,color = 'green')

#进行误差分析
SE = 0
for l in range(len(X)):
    SE = SE + (Y[[l]]-y[[l]])**2
    
RMSE = (float(SE)/(len(X)))**(0.5)

print("误差统计量SE=",SE)
print("       RMSE=",RMSE)
