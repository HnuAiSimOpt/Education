# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

# 二次方程
print("二次方程 y = a0 + a1x + a2x^2")
a0 = int(input("请输入系数: a0 = "))
a1 = int(input("请输入系数: a1 = "))
a2 = int(input("请输入系数: a2 = "))
# 随机点生成
X = np.arange(0, 5, 0.3)
Z = [ a0 + a1 * x + a2 * x ** 2 for x in X]
Y = np.array([np.random.normal(z,3) for z in Z])

plt.plot(X, Y, 'ro')
plt.show()

def get_data_a(X, Y):
    m = 3
    A = []
    # 每一个方程的系数
    for i in range(m):
        a = []
        # 当前方程中的每一个系数
        for j in range(m):
            a.append(sum(X ** (i+j)))
        A.append(a)
    return A

# 方程组的右端向量
def get_data_b(X, Y):
    m = 3
    b = []
    for i in range(m):
        b.append(sum(X**i * Y))
    return b

A = get_data_a(X, Y)
b = get_data_b(X, Y)

a0, a1, a2 = np.linalg.solve(A, b)

# 生成拟合曲线
_X = np.arange(0, 5, 0.1)
_Y = np.array([a0 + a1*x + a2*x**2 for x in _X])

plt.plot(X, Y, 'ro', _X, _Y, 'b', linewidth=2)
plt.title("y = {} + {}x + {}$x^2$ ".format(a0, a1, a2))
plt.show()

'''
# 线性方程
print("线性方程 y = a0 + a1x")
a0 = int(input("请输入系数: a0 = "))
a1 = int(input("请输入系数: a1 = "))

X = np.arange(0, 5, 0.3)
Z = [a0 + a1 * x for x in X]
Y = [np.random.normal(z, 0.5) for z in Z]

plt.plot(X, Y, 'ro')
plt.show()

def get_data(x, y):
    N = len(x)
    sumx = sum(x)
    sumy = sum(y)
    sumx2 = sum(x**2)
    sumxy = sum(x*y)
    A = np.mat([[N, sumx], [sumx, sumx2]])
    b = np.array([sumy, sumxy])
    return np.linalg.solve(A, b)

a0, a1 = get_data(X, Y)

_X = [0, 5]
_Y = [a0 + a1 * x for x in _X]

plt.plot(X, Y, 'ro', _X, _Y, 'b', linewidth=2)
plt.title("y = {} + {}x".format(a0, a1))
plt.show()
'''