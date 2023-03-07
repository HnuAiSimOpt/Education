import numpy as np
from scipy import linalg
from scipy.optimize import leastsq


def func(x, p):
    A, k, theta = p
    return A * np.sin(2 * np.pi * k * x + theta)


def residuals(p, y, x):
    return y - func(x, p)


x = np.linspace(0, -2 * np.pi, 100)
A, k, theta = 10, 0.34, np.pi / 6                     # 真实数据的函数参数
y0 = func(x, [A, k, theta])                           # 真实数据
y1 = y0 + 2 * np.random.randn(len(x))                 # 加入噪声之后的实验数据
p0 = [7, 0.2, 0]                                      # 第一次猜测的函数拟合参数

plsq = leastsq(residuals, p0, args=(y1, x))
print("真实参数:", [A, k, theta])
print("拟合参数:", plsq[0])

import matplotlib.pyplot as plt
import pylab as pl

plt.rcParams['font.sans-serif'] = ['SimHei']          # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False            # 用来正常显示负号
pl.plot(x, y0, marker='+', label=u"真实数据")
pl.plot(x, y1, marker='^', label=u"带噪声的实验数据")
pl.plot(x, func(x, plsq[0]), label=u"拟合数据")
pl.legend()
pl.show()
