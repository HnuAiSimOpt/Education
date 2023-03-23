import numpy as np
import matplotlib.pyplot as plt

# 计算最小二乘拟合直线的函数
def least_squares_fit(x, y):
    """
    :param x: 自变量 x 的值列表
    :param y: 因变量 y 的值列表
    :return: 拟合直线的截距和斜率
    """
    n = len(x)
    sum_x = sum(x)
    sum_y = sum(y)
    sum_x_square = sum([xi**2 for xi in x])
    sum_xy = sum([x[i]*y[i] for i in range(n)])
    # 计算矩阵的逆
    a = np.array([[n, sum_x], [sum_x, sum_x_square]])
    b = np.array([sum_y, sum_xy])
    a_inv = np.linalg.inv(a)
    # 计算截距和斜率
    c = np.dot(a_inv, b)
    intercept, slope = c[0], c[1]
    return intercept, slope

# 示例数据
x = [1, 2, 3, 4, 5]
y = [1.2, 1.9, 3.2, 4.1, 5.3]

# 调用函数计算拟合直线的截距和斜率
intercept, slope = least_squares_fit(x, y)

# 绘制散点图和拟合直线
plt.scatter(x, y, color='blue')
plt.plot(x, [slope * xi + intercept for xi in x], color='red')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Least squares fit line')
plt.show()

# 打印结果
print(f"拟合直线为 y = {slope:.2f}x + {intercept:.2f}")
