# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:06:09 2023
@author: 24211
"""

"""
此次最小二乘拟合为多项式拟合，项数为3，方法为矩阵求解法
首先引入所需库
"""
from matplotlib import pyplot as plt
import random
import numpy as np
"""
第一步，下方代码用于产生拟合样本数据和计算均方偏差数据
各20组"""
a0=random.uniform(1,5)
a1=random.uniform(1,5)
a2=random.uniform(1,5)
x=np.random.uniform(0,10,20)
e=np.random.normal(0,20,20)     #引入噪声
x1=np.random.uniform(0,10,20)   
e1=np.random.normal(0,20,20)
y=a0+a1*x+a2*x**2+e     #获得拟合数据样本
y1=a0+a1*x1+a2*x1**2+e1     #获得测验均方偏差样本
plt.scatter(x,y)  #在图中表达
np.set_printoptions(precision=2,suppress=True)
print('\n用于求拟合曲线的样本数据为\nx:{}\ny:{}'.format(x,y))
print('\n用于求均方偏差的样本数据为\nx:{}\ny:{}'.format(x1,y1))
"""
第二步，定义下方函数用于获得系数矩阵G和d，方便构造法方程，解法方程可得多项式系数
"""
def fai(j,k):
    s=0
    if j<=3:
        for i in range(0,10):
            faij=x[i]**j
            faik=x[i]**k
            s+=faij*faik
        return s    #此用于构造G
    else:
        for i in range(0,10):
            faik=x[i]**k
            s+=faik*y[i]
        return s    #此用于构造d
G=[[fai(0,0),fai(0,1),fai(0,2)],[fai(1,0),fai(1,1),fai(1,2)],[fai(2,0),fai(2,1),fai(2,2)]]
#获得系数矩阵G#
d=[[fai(10,0)],[fai(10,1)],[fai(10,2)]]
#获得矩阵d
G=np.array(G)
d=np.array(d)
a_=np.dot(np.linalg.inv(G),d)   #利用法方承Ga=d，得多项式系数向量a=G**-1 *d,从而解得拟合多项式系数向量
x_=np.arange(0,10,0.1)
s_=a_[0]+a_[1]*x_+a_[2]*x_**2   #得到拟合函数
plt.plot(x_,s_)#作拟合函数图
plt.show()
print('拟合函数\n第一项a0是{:.2f}\n第二项a1是{:.2f}\n第三项a2是{:.2f}'.format(float(a_[0]),float(a_[1]),float(a_[2])))
print('即s={:.2f}+{:.2f}x+{:.2f}x**2'.format(float(a_[0]),float(a_[1]),float(a_[2])))
"""
第三步 利用上述样本数据，求均方偏差MSE
"""
for i in range(20):
    MSE=0
    s=a_[0]+a_[1]*x1[i]+a_[2]*x1[i]**2
    p=(s-y1[i])**2
    MSE=MSE+p
print('均方偏差为{:.2f}'.format(float(MSE/20)))
    






