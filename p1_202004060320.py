#车辆2003  晏邦   202009060426   2022/3/8
#通过法方程和关系式Ga=d可以得到一个二元一次方程组，再通过克拉默公式得到拟合系数。WLS
import matplotlib.pyplot as plt
import numpy as np
import random 

#构造数据集
x=np.array([ 1.2,2.4 ,3.6 ,4.8 ,6.2 ,6.6,8.8,9.9,12,15.3])
y=1.48*x+2.56+random.random()       #假设实际的系数B0=1.48，B1=2.56

#构造拟合函数
def Linear_Fitting(x , y):      
    a00=1                         
    a01,a10,a11,b1,b2=0,0,0,0,0
    N = len(x)
    for k in range(0,N):                   
        a01=a10=a10+x[k]
    for m in range(0,N):
        a11=a11+x[m]**2
    for n in range(0,N):
        b1=b1+y[n]
    for i in range(0,N):
        b2=b2+x[i]*y[i]
    a0=a00*a11-a01*a10                      
    a1=b1*a11-a01*b2
    a2=a00*b2-a10*b1
    B0=-a1/a0
    B1=a2/a0
    return B0,B1
#样本检测
B0,B1=Linear_Fitting(x, y)         
y1=B0+B1*x
print("所求拟合曲线为P（x）=",B0,"+",B1,"x")
x1=np.array([2.8,5.9,8.5,10.4,13,16.4])            
y2=B1*x1+B0                              
y3=1.48*x+2.56
e=0.003*y                         
for i in range(0,len(x1)):
   while ((y2[i]-y3[i])**2) > e[i]:      
        B0 = B0-0.005
        B1 = B1-0.005
        y2 = B1*x1+B0
   while (abs(y2[i]-y3[i])-abs(y2[i]-y3[i])) <0:
       B0 = B0+0.005
       B1 = B1+0.005

plt.figure(figsize=(8, 6))
plt.scatter(x, y, color='red', label='Sample data', linewidth=2)   #绘制数据集的散点图
plt.plot(x, y1, color='blue', label='Fitting Curve', linewidth=2)  #绘制拟合后函数的直线图
plt.legend()  
plt.show()