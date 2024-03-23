# -*- coding: utf-8 -*-
"""
Spyder 编辑器

车辆2102 刘熠 202104060806
"""
import matplotlib.pyplot as plt
import numpy as np 
def Polynomial(k,targs):                 
    def arr_func(x):
        arr=1
        for i in targs:        
            if i==k: continue   
            arr*=x-i[0]
            arr/=k[0]-i[0]
        arr*=k[1]
        return arr
    return arr_func
def Lagrange(*targs):                  
    funcs=[Polynomial(i,targs) for i in targs]   
    def arr_func(x):
        arr=0
        for i in funcs:arr+=i(x)
        return arr
    return arr_func
data=[[1,1],[2,4],[3,9],[4,2],[5,6],[6,8]]
func=Lagrange(*data)           
x=np.arange(0,6,0.1)    
y=[func(i) for i in x]  
plt.title("Lagrange Interpolation Polynomial")
plt.plot(x,y)
plt.show()              