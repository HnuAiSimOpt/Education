'''蔺柯通
202078020226
2023.3.9'''

from numpy import *
import matplotlib.pyplot as plt
from pylab import mpl

"""一元线性拟合
采用的拟合数据为xi=1,2,3,4,5,6,7,8
对应的相应函数值yi=1.2,3.0,2,4,3.7,6,5.3,7.1
"""

x = [1, 2, 3, 4, 5, 6, 7, 8];
y = [1.2, 3.0, 2, 4, 3.7, 6, 5.3, 7.1]
x_ = mean(x)
y_ = mean(y)
l = len(x)
"""完成拟合曲线参数计算"""
i = 0
k1 = 0
k2 = 0
for i in range(l):
    k1 = k1 + x[i] * y[i]
    k2 = k2 + x[i] * x[i]
    i = i + 1
k = (k1 - l * x_ * y_) / (k2 - l * x_ * x_)
b = y_ - k * x_

"""完成完后曲线上相应的函数值的计算"""
i = 0
y1 = []
for i in range(l):
    y1.append(k * x[i] + b)
    i = i + 1

"""完成函数的绘制"""


def draw(data_x, data_y_new, data_y_old):
    plt.plot(x, y1, label="拟合曲线", color="black")
    plt.scatter(x, y, label="离散数据")
    mpl.rcParams['font.sans-serif'] = ['SimHei']
    mpl.rcParams['axes.unicode_minus'] = False
    plt.title("一元线性拟合数据")
    plt.legend(loc="upper left")
    plt.show()


draw(x, y1, y)

'''进行误差分析'''
SE = 0
i = 0
for i in range(l):
    SE = SE + (y1[i] - y[i]) ** 2
    i = i + 1

print("误差统计量SE=", SE)