# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 16:31:59 2023

@author: 王熙华 202004060832 车辆2004班
"""
#date：2023.3.6
#purposal:实际函数为y=x**2+1,用最小二乘法对其进行一次线性拟合

import random

size1 = 100 #进行线性拟合点的数量
size2 = 1000 #进行误差分析的点的数量
sumx = 0 #100个点的x之和
sumy = 0 #100个点的y之和
sumxy = 0 #100个点的x*y之和
sumx2 = 0 #100个点的x**2之和
#建立一个列表ALL存储1200个随机数
ALL = []
for i in range(1200):
    a1 =random.random()*10
    ALL.append(a1)
#从ALL列表中取出1100个不同的随机数并将前100 个放入列表X1，后1000个放入列表X2
ALLX = random.sample(ALL,1100)
X1 = ALLX[0:100]
X2 = ALLX[101:]
#通过循环来求解一次拟合函数的a，b 并输出拟合函数
for i in X1:
    sumx += i
    sumy += (i*i+1)
    sumxy += (i*i*i+i)
    sumx2 += (i*i)
avaragex = sumx/size1
avaragey = sumy/size1
b = (sumxy - size1 * avaragex *avaragey)/(sumx2 - size1 * avaragex *avaragex)
a = avaragey - b * avaragex
print("y = {0}*x + ({1})".format(b, a))
#代入X2中的值对拟合函数和原函数进行误差分析并输出
for x in X2:
    y = b * x + a
    error = abs(y) - x*x -1
    print("%f,%f,%f"%(x,y,error))



   

    
    



