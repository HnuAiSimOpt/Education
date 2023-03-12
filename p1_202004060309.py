# -*- coding: utf-8 -*-
"""
author: 朱晨侃 202004060309 车辆2004班
date：2023.3.9
purposal:OLS拟合
"""
import random

sz1 = 100 #进行线性拟合点的数量
sz2 = 1000 #进行误差分析的点的数量
sumx = 0 #定义点的 x 之和初值为0
sumy = 0 #定义点的 y 之和初值为0
sumxy = 0 #定义点的 x*y 之和初值为0
sumx2 = 0 #定义点的 x**2 之和初值为0
#建立一个列表AX存储1100个随机数作为随机点的x值
AX = []
for i in range(1100):
    a1 =random.random()*100
    AX.append(a1)
#从AX列表中取出1100个不同的随机数并将前100 个放入列表X1，后1000个放入列表X2
X1 = AX[0:100]
X2 = AX[100:]
#建立一个列表AY存储1100个随机数作为随机点的y值
AY = []
for i in range(1100):
    a1 =random.random()*100
    AY.append(a1)
#从AY列表中取出1100个不同的随机数并将前100 个放入列表Y1，后1000个放入列表Y2
Y1 = AY[0:100]
Y2 = AY[100:]
#通过循环来求解一次拟合函数的a，b，并输出拟合函数
for i in range(0,99):
    sumx += X1[i]
    sumy += Y1[i]
    sumxy += X1[i] * Y1[i]
    sumx2 += X1[i] * X1[i]
avaragex = sumx/sz1
avaragey = sumy/sz1
b = (sumxy - sz1 * avaragex *avaragey)/(sumx2 - sz1 * avaragex *avaragex)
a = avaragey - b * avaragex
print("y = {0}*x + ({1})".format(b, a))
#用拟合函数求出y值并与实际点进行均方根误差计算，并输出误差值
sumd2 = 0
for i in range(0,999):
    y = b * X2[i] + a
    sumd2 += (y-Y2[i])*(y-Y2[i])
error = (sumd2 / 1000)**0.5
print(error)