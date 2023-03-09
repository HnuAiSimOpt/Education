# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import math
import random

#构造拟合用函数
def linefit(x , y):
  N = float(len(x))
  sumx,sumy,sumxx,sumyy,sumxy=0,0,0,0,0
  for i in range(0,int(N)):
    sumx += x[i]
    sumy += y[i]
    sumxx += x[i]*x[i]
    sumyy += y[i]*y[i]
    sumxy += x[i]*y[i]
  a = (sumy*sumx/N -sumxy)/( sumx*sumx/N -sumxx)
  b = (sumy - a*sumx)/N
  r = abs(sumy*sumx/N-sumxy)/math.sqrt((sumxx-sumx*sumx/N)*(sumyy-sumy*sumy/N))
  return a,b,r

#数据样本
x=np.array([ 1.1 ,2.2 ,3.3 ,4.4 ,5.5 ,6.6])
y=2*x+1.5+random.random()
a,b,r=linefit(x,y)                  #计算出a，b，r的值
y1=a*x+b                            #拟合的直线
x1=np.array([6,7,8,9,10,11])        #测试点
y2=a*x1+b                           #测试点模拟函数值
y3=2*x1+1.5                         #测试点真实值
e=0.01*y3                           #给予参考误差值

#绘制图表
plt.figure(figsize=(8, 6))
plt.scatter(x, y, color='red', label='Sample data', linewidth=2)
plt.plot(x, y1, color='blue', label='Fitting Curve', linewidth=2)
plt.scatter(x1, y3, color='yellow', label='Sample data', linewidth=2)
plt.scatter(x1, y2, color='green', label='Sample data', linewidth=2)
plt.plot(x1, y2, color='purple', label='Fitting Curve', linewidth=2)
plt.legend()
plt.xlabel('x', fontproperties='simHei', fontsize=12)
plt.ylabel('y', fontproperties='simHei', fontsize=12)
plt.show()

print(a,b)