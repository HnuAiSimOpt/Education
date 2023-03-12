# coding=utf-8  
''''' 
车辆2003 
蔡江羽
2023/3/8  
202004061332
'''  
import matplotlib.pyplot as plt    
import numpy  
import random  
  
fig = plt.figure()  
ax = fig.add_subplot(111)  #绘制画布
  
order=1  #阶数为1阶  
  
#生成曲线上的各个点  
x = numpy.arange(-1,1,0.1)  #生成x轴数据
y = [5*a+6 for a in x]  #生成y轴数据

#偏移各点，并放入到xa,ya中去  
i=0  
xa=[]  
ya=[]  
for xx in x:  
    yy=y[i]  
    d=float(random.randint(80,120))/100   
    i+=1  
    xa.append(xx*d)  
    ya.append(yy*d)  
ax.plot(xa,ya,color='r',linestyle='',marker='.')  
  
#进行曲线拟合  
matA=[]  
for i in range(0,order+1):  
    matA1=[]  
    for j in range(0,order+1):  
        tx=0.0  
        for k in range(0,len(xa)):  
            dx=1.0  
            for l in range(0,j+i):  
                dx=dx*xa[k]  
            tx+=dx  
        matA1.append(tx)  
    matA.append(matA1)  
matA=numpy.array(matA)  
  
matB=[]  
for i in range(0,order+1):  
    ty=0.0  
    for j in range(0,len(xa)):  
        dy=1.0  
        for k in range(0,i):  
            dy=dy*xa[j]  
        ty+=ya[j]*dy  
    matB.append(ty)  
   
matB=numpy.array(matB)  
  
matAA=numpy.linalg.solve(matA,matB)  
  
#画出拟合后的曲线  
xxa= numpy.arange(-1,1,0.01)  
yya=[]  
for i in range(0,len(xxa)):  
    yy=0.0  
    for j in range(0,order+1):  
        dy=1.0  
        for k in range(0,j):  
            dy*=xxa[i]  
        dy*=matAA[j]  
        yy+=dy  
    yya.append(yy)  
ax.plot(xxa,yya,color='g')  
  
plt.show() 
