# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 13:31:18 2024

@author: 25600
"""
#  陈毅 202104060829 车辆2102
# 定义拉格朗日函数
def lagrange(x,y,Xn):                           #x为已知插值点数值
    Ln=0                                        #y为已知插值点对应函数值
    for i in range(len(x)):                     #Xn为预计算插值点
        L=y[i]
        for j in range(len(x)):
            if i!=j:
                L*=(Xn-x[j])/(x[i]-x[j])
        Ln+=L
    return round(Ln,3)                          #保留3位小数
x=eval(input("请输入 x="))
y=eval(input("请输入 y="))
Xn=eval(input("请输入 Xn="))
Ln=lagrange(x,y,Xn)
print("Ln=",Ln)


# 例子
#请输入 x=1,3,5
#请输入 y=2,10,1
#请输入 Xn=10
#得出新插值点为Ln=-95.875

