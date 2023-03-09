
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from numpy.linalg import solve
from itertools import chain
n = 20                      #总共20个数
x = np.linspace(0, 20, n)   #20个数在0-20的区间均匀分布
a = np.random.normal(size=n)#利用random随机生成偏差a
y = x + a                  #y等于x值加上偏差
sum1=sum(x)                #xi的和
sum2=sum(y)                #yi的和
sum3=sum(np.multiply(np.array(x),np.array(x)))    #xi平方的和
sum4=sum(np.multiply(np.array(x),np.array(y)))    #xi与yi积的和
H = np.array([[10,sum1], [sum1,sum3]])          #产生格拉姆矩阵
h = np.array([sum2,sum4])
solution =solve(H,h)                                #用solve得回归方程的解
y_fit=solution[0]+solution[1]*x                     #回归方程
print('回归方程：y={0}+{1}x'.format(solution[0],solution[1]))#写出回归方程
fig, ax = plt.subplots(figsize=(18,14))             #生成Figure画布和axes对象
ax.plot(x, y, 'o')                                  #用ax对象在区域内生成点图
ax.plot(x, y_fit, 'r--.')                           #用ax对象在区域内生成曲线
plt.show()