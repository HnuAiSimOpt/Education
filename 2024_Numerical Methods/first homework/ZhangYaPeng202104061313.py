# 张亚鹏 202104061313
import matplotlib.pyplot as plt
import numpy as np

# 定义拉格朗日函数
def lagrange(x1, y1, x_values):
    n = len(x1)
    y_values = []
    for x in x_values:
        y = 0
        for k in range(n):
            t = 1
            for i in range(n):
                if i != k:
                    t *= (x - x1[i]) / (x1[k] - x1[i])
            y += t * y1[k]
        y_values.append(y)
    return y_values

# 插值节点
x1 = [0, 1, 2, 3, 4, 5]
y1 = [0, 2.2, 4.0, 5.8, 7.1, 7.3]

# 插值点
xx = [2.5, 3.5]

# 调用拉格朗日插值函数,输出插值点xx的函数值yy
yy = lagrange(x1, y1, xx)
for i, j in zip(xx, yy):
    print(f"插值点 x = {i} 的函数值 y = {j}")

# 绘制图像
plt.scatter(x1, y1, color='red', label='Interpolation nodes')
plt.scatter(xx, yy, color='blue', label='Output points')
x_continuous = np.linspace(min(x1), max(x1), 100)
y_continuous = lagrange(x1, y1, x_continuous)
plt.plot(x_continuous, y_continuous, label='Lagrange Interpolation')
plt.legend()
plt.show()

