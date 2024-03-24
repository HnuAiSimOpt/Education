#赵天乐 202109060504

import numpy as np  
  
def lagrange(datax,datay,x):  
    n = len(datax)  
    L = np.zeros(n)   
    y = 0  
    for i in range(n):  
        L[i] = 1  
        for j in range(n):  
            if i != j:  
                L[i] *= (x-datax[j])/(datax[i]-datax[j])  
        y += L[i]*datay[i]  
  
    return y  
datax= np.array(input("请输入数据点的x坐标，用空格分隔：").split(), dtype=float)  
datay= np.array(input("请输入数据点的y坐标，用空格分隔：").split(), dtype=float)  
x = float(input("请输入需要插值的x坐标："))  
y = lagrange(datax,datay, x)  
print(f"插值结果: y = {y}")