# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 16:15:17 2023

@author: your daddy
"""
import numpy as np
import math
import random 
#首先做出拟合函数
def niheline(x,y):
  n=float(len(x))
  sum_x=0 #初始化
  sum_y=0
  sumx_x=0
  sumy_y=0
  sum_xy=0
  for i in range(0,int(n)):
    sum_x+=x[i]#求x和
    sum_y+=y[i]#求y和
    sumx_x+=x[i]*x[i]#求x平方和
    sumy_y+=y[i]*y[i]#求y平方和
    sum_xy+=x[i]*y[i]#求xy和
  a=(sum_y*sum_x/n-sum_xy)/(sum_x*sum_x/n-sumx_x)
  b=(sum_y-a*sum_x)/n
  er=abs(sum_y*sum_x/n-sum_xy)/math.sqrt((sumx_x-sum_x*sum_x/n)*(sumy_y-sum_y*sum_y/n))
  return a,b,er
#设点并验证函数拟合
x=np.array([2,3,4,5,6,7,8,11]) #给出样本点
y=3*x+4+random.uniform(0,2) #拟合曲线方程
a,b,er=niheline(x,y) #算出拟合曲线方程参数
y1=a*x+b
x1=np.array([10,13,16,19,22,25,28,37]) #给出精准数据
y2=a*x1+b
y3=3*x1+4 #原曲线方程
standard=y3*0.02 #误差设定
for i in range(0,len(x1)): #依据最小二乘法降低误差
   while((y2[i]-y3[i])**2)>standard[i]:
        a=a-0.01
        b=b-0.01
        y2=a*x1+b
#输出拟合函数斜率和截距
print(a,b)