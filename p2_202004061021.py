"""本文python代码实现的是最小二乘法线性拟合
问题：对直线附近的带有噪声的数据进行线性拟合，最终求出w,b的估计值。
最小二乘法基本思想是使得样本方差最小。
代码中self_func()函数为自定义拟合函数"""

"""
此次最小二乘拟合为线性拟合，首先引入所需库
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

"""
第一步，下方代码用于产生拟合样本数据和计算均方偏差数据
各100组"""

n = 101
 
x = np.linspace(0,10,n)
noise = np.random.randn(n)
y = 4.5 * x + 0.6 + 2.0 * noise

"""
第二步，定义自定义拟合函数self_func()，用于对样本数据进行线性拟合，求出w,b的估计值，计算均方偏差数据。
"""
 
def self_func(steps=100, alpha=0.01):
  w = 0.5
  b = 0
  alpha = 0.01
  for i in range(steps):
    y_hat = w*x + b
    dy = 2.0*(y_hat - y)
    dw = dy*x
    db = dy
    w = w - alpha*np.sum(dw)/n
    b = b - alpha*np.sum(db)/n
    e = np.sum((y_hat-y)**2)/n
  print ('self_func:\tW =',w,'\n\tb =',b,'\n\te =',e)
  plt.scatter(x,y)
  plt.plot(np.arange(0,10,1), w*np.arange(0,10,1) + b, color = 'r', marker = 'o', label = 'self_func(steps='+str(steps)+', alpha='+str(alpha)+')')
 
"""
第三步 利用上述定义函数，输出w,b的值和均方偏差e
"""
self_func(10000)
plt.legend(loc='upper left')
plt.show()
 