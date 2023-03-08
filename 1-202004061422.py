#车辆2003 杨明焕 202004061422  2023/3/7
#通过输入一系列x的值，然后在理论函数下得到一系列的点。然后通过法方程和关系式Ga=d可以得到一个二元一次方程组，再通过克拉默公式得到拟合系数。
#这种方法的缺陷是需要解方程组，按照我现在的知识水平只会用克拉默法则，因此能拟合的阶数不高。相信通过后面迭代课程的学习，我可以进行更高阶数的拟合

import matplotlib.pyplot as plt
import numpy as np                      #这里导入numpy库只是为了画图
#第一步：通过一个理论函数建立一个数据集
reality_x = []
reality_y = []
print("please input sample's number")    #输入样本的个数
N=int(input())
print('please input sample x')           #输入样本
for i in range(0,N):
    reality_x.append(float(input()))
B0=1.5                                   #假设实际的系数B0=1.5，B1=1.2
B1=1.2
for j in range(0,N):                     #根据输入的x和假设的函数得到数据集
    reality_y.append(B0+B1*reality_x[j])
      
#第二步：通过法方程公式和克拉默法则(因为只有一个2*2的矩阵，比较简单）得到拟合的系数Fit_B0和Fit_B1
a00=1
a01=0
a10=0
a11=0
b1=0
b2=0
for k in range(0,N):                   #法方程的计算，法方程是多元函数对每一个系数求偏导数以后得到的线性方程组
    a01=a10=a10+reality_x[k]
for m in range(0,N):
    a11=a11+reality_x[m]**2
for n in range(0,N):
    b1=b1+reality_y[n]
for l in range(0,N):
    b2=b2+reality_x[l]*reality_y[l]
a0=a00*a11-a01*a10                     #用克拉默法则求解线性方程组
a1=b1*a11-a01*b2
a2=a00*b2-a10*b1
Fit_B0=-a0/a1
Fit_B1=a0/a2
print("所求拟合曲线为P（x）=",Fit_B0,"+",Fit_B1,"x")

#第三步：样本点检验和误差分析
print('please input text sample x')           #输测试入样本
text_x=float(input())
e=B0+B1*text_x-Fit_B0-Fit_B1*text_x
print("输入的测试点误差为：")
print(e)

#第四步：绘图
plt.figure()
plt.scatter(reality_x, reality_y, c="red",label="reality")
plt.xlabel("x")
plt.ylabel("y")
x=np.linspace(0,100)
y=B0+B1*x
plt.plot(x,y,c="b")
plt.show()



