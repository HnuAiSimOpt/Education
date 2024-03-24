# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 20:57:02 2024

@author: 王璇 202104060120
"""
import numpy as np
#键入节点和函数值
X0=np.array(list(map(int,input("请输入插值节点，用空格分割：").split())))
Y0=np.array(list(map(int,input("请输入函数值，用空格分割：").split())))
#节点数
n=len(X0)
#拟合坐标点横坐标数组，使用数组同时计算拟合曲线上一千个点的坐标
x=np.linspace(0,X0[-1],1000)
#初始化插值多项式
Ln=np.zeros(len(x))
#使用sympy得符号表达式
import sympy
x1=sympy.Symbol('x')
Lnbd=0
#循环计算不同的插值基函数，同时得到符号表达式
for i in range(n):
#初始化插值基函数
    lk=np.ones(len(x))
    lkbd=1
#根据δ条件计算插值基函数
    for j in range(n):
        if i!=j:
            lk*=(x-X0[j])/(X0[i]-X0[j])
            lkbd*=(x1-X0[j])/(X0[i]-X0[j])
#插值基函数线性组合成插值多项式
    Ln+=Y0[i]*lk
    Lnbd+=Y0[i]*lkbd
print("拉格朗日插值多项式为Ln={}".format(Lnbd))
# 绘图
import matplotlib.pyplot as plt
plt.plot(X0,Y0,'ro',label='row')
plt.plot(x,Ln,'b-',label='lagrange')
plt.legend()
plt.show()
