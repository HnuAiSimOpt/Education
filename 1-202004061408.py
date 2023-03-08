import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from numpy.linalg import solve
from itertools import chain
n = 10
x = np.linspace(0, 10, n)   #生成10个数在0-10之间均匀排列
e = np.random.normal(size=n)#按正态分布生成偏差
y = x + e                   #y的实际值
                            #以下是书上例9求解拟合曲线的方法
sum_x=sum(x)                #所有x的和
sum_y=sum(y)                #所有y的和
sum_xx=sum(np.multiply(np.array(x),np.array(x)))    #所有xi平方的和
sum_xy=sum(np.multiply(np.array(x),np.array(y)))    #所有xi乘yi的和
G = np.array([[10,sum_x], [sum_x,sum_xx]])          #求解方程组，生成格拉姆矩阵
d = np.array([sum_y,sum_xy]) 
solution =solve(G,d)                                #方程的解
y_fit=solution[0]+solution[1]*x                     #回归方程
print('回归方程为：y={0}+{1}x'.format(solution[0],solution[1]))
fig, ax = plt.subplots(figsize=(16,12))             #绘制数据图和回归方程曲线
ax.plot(x, y, 'o', label='数据')
ax.plot(x, y_fit, 'r--.',label='拟合结果')
plt.show()
