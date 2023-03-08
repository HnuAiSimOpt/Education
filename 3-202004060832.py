# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 16:31:59 2023

@author: 王熙华 202004060832 车辆2004班
"""
#date：2023.3.6
#purposal:找100个随机点,用最小二乘法对其进行一次线性拟合，然后用1000个样本点进行均方根误差分析

import random

size1 = 100 #进行线性拟合点的数量
size2 = 1000 #进行误差分析的点的数量
sumx = 0 #100个点的x之和
sumy = 0 #100个点的y之和
sumxy = 0 #100个点的x*y之和
sumx2 = 0 #100个点的x**2之和
#建立一个列表ALLX存储1200个随机数作为随机点的x值
ALLX = []
for i in range(1200):
    a1 =random.random()*10
    ALLX.append(a1)
#从ALLX列表中取出1100个不同的随机数并将前100 个放入列表X1，后1000个放入列表X2
ALLx = random.sample(ALLX,1100)
X1 = ALLx[0:100]
print(len(X1))
X2 = ALLx[100:]
print(len(X2))
#建立一个列表ALLY存储1200个随机数作为随机点的y值
ALLY = []
for i in range(1200):
    a1 =random.random()*10
    ALLY.append(a1)
#从ALLY列表中取出1100个不同的随机数并将前100 个放入列表Y1，后1000个放入列表Y2
ALLy = random.sample(ALLY,1100)
Y1 = ALLy[0:100]
Y2 = ALLy[100:]
#通过循环来求解一次拟合函数的a，b 并输出拟合函数
for i in range(100):
    sumx += X1[i]
    sumy += Y1[i]
    sumxy += X1[i] * Y1[i]
    sumx2 += X1[i] * X1[i]
avaragex = sumx/size1
avaragey = sumy/size1
b = (sumxy - size1 * avaragex *avaragey)/(sumx2 - size1 * avaragex *avaragex)
a = avaragey - b * avaragex
print("拟合函数为：y = {0}*x + ({1})".format(b, a))
#用拟合函数求出y值并与实际点进行均方根误差分析，并输出误差值
sumd2 = 0
for i in range(1000):
    y = b * X2[i] + a
    sumd2 += (y-Y2[i])*(y-Y2[i])
error = (sumd2 / 1000)**0.5
print("拟合函数均方根误差为：",format(error))
    
   



   

    
    



