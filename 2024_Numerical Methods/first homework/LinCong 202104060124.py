#车辆2101 202104060124 林聪
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange

# 定义原多项式
def f(x):
    return x**4*7-4*x**3+7

# 拉格朗日插值函数
def lagrange_interpolation(x, y, val):
    n = len(x)
    y_interp = 0
    for i in range(n):
        L = 1
        for j in range(n):
            if j != i:
                L *= (val - x[j]) / (x[i] - x[j])
        y_interp += y[i] * L
    return y_interp

# 已知数据点
x_points = np.array([-100, -50, 5, 60, 100])
y_points = f(x_points)

# 插值点
x_interp = np.linspace(min(x_points), max(x_points), 201)
y_interp = np.array([lagrange_interpolation(x_points, y_points, xi) for xi in x_interp])

'''库验证
print(lagrange(x_points, y_points))
print(lagrange(x_points, y_points)(x_interp))''' 

# 绘制原多项式
plt.plot(x_interp, f(x_interp), label='Original Polynomial', color='blue')

# 绘制插值多项式
plt.plot(x_interp, y_interp, label='Lagrange Interpolation', color='red', linestyle='--')

# 设置图例
plt.legend()

# 显示图形
plt.show()


