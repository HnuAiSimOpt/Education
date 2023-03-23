#用最小二乘法对y=x**(0.5)+x进行一次线性拟合
import random
L1 = []
for i in range(500):
    m =random.random()*10
    L1.append(m)            #建立一个列表L1存储500个用于线性拟合的随机数
L2 = []
for i in range(1000):
    m =random.random()*10
    L2.append(m)            #建立一个列表L2存储1000个用于线性拟合的随机数
sumx=sum(L1)                #500个点的x之和
sumy = 0;sumxy = 0;sumx2 = 0#500个点的y/x*y/x**2之和
for i in L1:
    sumy += (i**(0.5)+i)
    sumxy += (i**(0.5)*i+i*i)
    sumx2 += (i*i)
avex=sumx/len(L1)           #x的平均值
avey=sumy/len(L1)           #y的平均值
k= (sumxy - len(L1) * avex *avey)/(sumx2 - len(L1) * avex **2)
b= avey-k*avex
print("y = {0}*x + ({1})".format(k, b))
#误差分析
L3=[]
for x in L2:
    y = k * x + b
    error = y - x**(0.5) -x
    L3.append(error)
    print("%f,%f,%f"%(x,y,error))
print(sum(L3))
