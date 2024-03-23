# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 13:52:15 2024

@author: akai
"""
#朋锴 车辆2102 202104060707
import numpy as np 
# 定义拉格朗日插值法函数
def lagrange_interpolation(x, y, x_xin):
    # x: 已知数据点
    # y: 已知数据点对应的函数值
    # x_xin: 需要插值的新数据点
    
    # 定义插值后的函数值
    y_xin = 0
    
    # 查找每种数据点的个数并进行循环，计算拉格朗日多项式
    for i in range(len(x)):
        # 初始化基函数为1
        L = 1
        
        # 计算拉格朗日基函数
        for j in range(len(x)):
            if i != j:
                L *= (x_xin - x[j]) / (x[i] - x[j])
        
        # 加权求和得到拉格朗日插值多项式的值
        y_xin += y[i] * L
    
    return y_xin
 
# 例如求出如下插值
x = np.array([1, 2, 3, 4, 5, 6])#插值节点
y = np.array([1, 3, 7, 15, 31, 63])#插值节点
x_xin = np.array([2.5, 3.5, 4.5, 5.5])#插值点
 
# 用拉格朗日插值法函数计算新点的函数值
for xin_x in x_xin:
    y_xin = lagrange_interpolation(x, y, xin_x)
    print(f"x = {xin_x}, y = {y_xin}")