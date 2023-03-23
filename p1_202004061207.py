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
  sum_x=0
  sum_y=0
  sumx_x=0
  sumy_y=0
  sum_xy=0
  for i in range(0,int(n)):
    sum_x+=x[i]
    sum_y+=y[i]
    sumx_x+=x[i]*x[i]
    sumy_y+=y[i]*y[i]
    sum_xy+=x[i]*y[i]
  a=(sum_y*sum_x/n-sum_xy)/(sum_x*sum_x/n-sumx_x)
  b=(sum_y-a*sum_x)/n
  er=abs(sum_y*sum_x/n-sum_xy)/math.sqrt((sumx_x-sum_x*sum_x/n)*(sumy_y-sum_y*sum_y/n))
  return a,b,er
#设点并验证函数拟合
x=np.array([2,3,4,5,6,7,8,11])
y=3*x+4+random.uniform(0,2)
a,b,er=niheline(x,y)
y1=a*x+b
x1=np.array([10,13,16,19,22,25,28,37])
y2=a*x1+b
y3=3*x1+4
standard=y3*0.02
for i in range(0,len(x1)):
   while((y2[i]-y3[i])**2)>standard[i]:
        a=a-0.01
        b=b-0.01
        y2=a*x1+b
#输出拟合函数斜率和截距
print(a,b)