import random
import matplotlib.pyplot as plt
import numpy as np
x = random.sample(range(0,200),60)#生成随机的60个真实点
x.sort()#对点进行排序
bo = 0
ao = 0
y = []
while True:##自定义设置目标实际函数
    if bo==0:
        bo=float(input('输入“b*x^2 + a”中的b，例如：4\n'))
    if ao==0:
        ao=float(input('输入“b*x^2 + a”中的a，例如：3\n'))
    if bo!=0 and ao!=0:
        break
for i in x:##取得实际函数对应的真实点
    j = bo*i**2 + ao
    y.append(j)
def gauss_noisy(x, y):##引用高斯噪声函数
    mu = 0
    sigma = 4
    for i in range(len(x)):
        x[i] += random.gauss(mu, sigma)
        y[i] += random.gauss(mu, sigma)
gauss_noisy(x, y)# 加入高斯噪声

###对函数进行拟合
##求解所需系数
len_x = len(x)
sum_x = np.sum(x)
sum_y = np.sum(y)
sum_xx = 0
sum_xy = 0
for i in range(len(x)):
    sum_xy = sum_xy + x[i]*y[i]
    sum_xx = sum_xx + x[i]**2
##开始求解a和b的值
A = np.array([[len_x,sum_x,sum_y],[sum_x,sum_xx,sum_xy]])
c = -A[1][0]/A[0][0]
A[1] = A[0]*c + A[1]
b = A[1][2]/A[1][1]
a = (A[0][2]-b*A[0][1])/A[0][0]
print("拟合直线a值为{:.2f}".format(a))
print("拟合直线b值为{:.2f}".format(b))

###误差估计（和反差sse）
xw = random.sample(range(0,200),100)##取100个点进行误差估计
xw.sort()
yw = []
for i in xw:
    j = bo*i**2 + ao
    yw.append(j)
gauss_noisy(xw, yw)
##另取得100个实际值
##求sse
sse = 0
for i in range(len(xw)):
    sj  = b * xw[i] + a
    sse = sse + (sj-yw[i])**2
print("估计的和方差为{:.2f}".format(sse))

##画出拟合曲线和离散点
X = np.linspace(0, 200, 1000)
Y = b * X + a
plt.plot(X, Y)
plt.title('Y={}+{}X'.format(a, b))
plt.legend(loc='lower right')
plt.scatter(x, y)# 画出这些加入噪声后的真实点
plt.show()
