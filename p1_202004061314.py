'''最小二乘法OLS车辆2003 202004061314 郝丁亮'''

import numpy as np
import matplotlib.pyplot as plt
X = np.arange(0, 15, 0.5)        #生成随机点
Z = [4 + 3 * x for x in X]      #原函数
Y = [np.random.normal(z, 0.5) for z in Z]

#计算x和y的平均值
n = len(X)
x_mean = sum((x) for x in X) / n
y_mean = sum((y) for y in Y) / n

#计算x和y的方差及协方差
var_x = sum((x - x_mean)**2 for x in X) / (n-1)
var_y = sum((y - y_mean)**2 for y in Y) / (n-1)
cov_xy = sum((x - x_mean)*(y - y_mean) for x,y in zip(X,Y)) / (n-1)

#求解
k = cov_xy / var_x
b = y_mean - k * x_mean

print("斜率:", k)
print("截距:", b)
_Y = [b+ k * x for x in X]
#红色圆点为原始数据，蓝色拟合直线
plt.plot(X, Y, 'ro', X, _Y, 'b', linewidth=2)
plt.title("y = {} + {}x".format(b, k))
plt.show()