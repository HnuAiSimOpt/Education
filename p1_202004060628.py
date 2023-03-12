#最小二乘法拟合散点图
#车辆2004邓婕 202004060628


import matplotlib.pyplot as plt
import numpy as np
plt.figure()

#设置初始数据
x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
y = [2, 4.1, 6.2, 7.8, 10, 11, 15, 16.2, 18.3, 19];

# 绘制散点
plt.scatter(x[:], y[:], 3, "red")

#计算斜率a与截距b
sum_x = 0;  # x的和
sum_y = 0;  # y的和
x_mean = 0  # x的均值
y_mean = 0  # y的均值
numerator = 0.0  # 分子
denominator = 0.0  # 分母
i = 0
while i < 10:
    sum_x += x[i];
    sum_y += y[i];
    i += 1;
x_mean = sum_x / 10;
y_mean = sum_y / 10;
for x_i, y_i in zip(x, y):
    numerator += (x_i - x_mean) * (y_i - y_mean)
    denominator += (x_i - x_mean) ** 2
a = numerator / denominator  # 求得a，斜率
b = y_mean - a * x_mean  # 求得b，截距
print('y=', '%.2f' % a, '*x+', '%.2f' % b)

#计算拟合后对应的函数值
x_test = np.linspace(1, 11, 10)
y_predict = a * x_test + b

# 计算和方差SSE
SSE = 0
j = 0
while j < 10:
    test_error = y_predict[j] - y[j]
    SSE += test_error ** 2
    j += 1
print("拟合函数的和方差为：", SSE)

#绘制拟合直线
plt.scatter(x, y, color='b')
plt.plot(x_test, y_predict, color='r')
plt.xlabel('x', fontproperties='simHei', fontsize=15)
plt.ylabel('y', fontproperties='simHei', fontsize=15)
plt.show()