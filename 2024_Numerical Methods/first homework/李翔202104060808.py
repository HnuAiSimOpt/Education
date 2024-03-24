# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 13:19:12 2024

@author: Lenovo
"""
#李翔 202104060808


x1=int(input('输入测试值x='))         #x为测试数值
n=int(input('输入插值点个数n='))          #n是插值点个数
XL=[]                   #XL为输入的x的插值
YL=[]                   #YL为输入的y的插值
for i in range(n):
    x=int(input('依次输入n个x的值='))
    XL.append(x)
for i in range(n):
    y=int(input('依次输入n个y的值='))
    YL.append(y)
def lagrange(x1):
    P=[]
    L=0
    for i in range(len(XL)):
        a=1
        b=1
        for j in range(len(XL)):
            if j!=i:
                a*=(x1-XL[j])
                b*=(XL[i]-XL[j])
        P.append(a/b)
    for i in range(len(YL)):
        L+=YL[i]*P[i]
    return L
print('结果为',lagrange(x1))

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    