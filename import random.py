import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

plt.rcParams['font.sans-serif'] = ['SimHei']
# 建立理论模型 y = x**2 在(0,100)内的随机100个样本库
Xi = np.array(random.sample(range(0, 100), 100))
Yi = np.array(Xi * Xi)  # 得到y = x方的对应样本数组


# 定义拟合函数的形式
def func(p, x):
    k, b = p
    return k * x + b


# 定义误差函数
def error(p, x, y, s):
    print(s)
    return func(p, x) - y


# 拟定参数的估计值
p0 = [100, 0]
# 用leastsq()函数进行参数估计
s = '参数估计次数'
Para = leastsq(error, p0, args=(Xi, Yi, s))
k, b = Para[0]
# 输出拟合函数及其参数
print("k=", k, "b=", b)
print("cost：" + str(Para[1]))
print("求解的拟合直线为:")
print("y=" + str(round(k, 2)) + "x+" + str(round(b, 2)))
# 求该线性函数相对于真值的均方误差
# 范围内 生成1000个样本点
X1 = np.array(np.linspace(0, 100, 1000))
Y1 = np.array(X1 * X1)
# 得到拟合函数的样本函数值
Y1_fit = np.array(k * X1 + b)
test_error = Y1 - Y1_fit
test_array = test_error ** 2
# 算出均方误差MSE
MSE = sum(test_array) / len(test_array)
print("该拟合函数的均方误差为:", round(MSE, 2))
# 图形可视化
plt.figure(figsize=(8, 6))
# 绘制参考数据点的散点图
# 作出样本点
plt.scatter(Xi, Yi, color="green", label="样本数据", linewidth=2)
plt.xlabel('x')
plt.ylabel('y')
x = np.linspace(0, 100, 1000)
y = k * x + b
# 做出拟合的直线
plt.plot(x, y, color="red", label="拟合直线", linewidth=2)
plt.title('y={}+{}x'.format(b, k))
plt.legend(loc='lower right')
plt.show()