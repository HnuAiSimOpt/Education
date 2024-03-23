#车辆2101班 桂楷森 202104060321
import numpy as np
import matplotlib.pyplot as plt

def input_points():
    x_values = input("请输入x值（用逗号分隔）: ").split(',')
    y_values = input("请输入对应的y值（用逗号分隔）: ").split(',')
    
    x_values = [float(x.strip()) for x in x_values]
    y_values = [float(y.strip()) for y in y_values]
    
    return x_values, y_values

def lagrange(x, x1, y1):
    L = 0
    for k in range(len(x1)):
        w = 1
        for m in range(len(x1)):
            if m != k:
                w *= (x - x1[m]) / (x1[k] - x1[m])
        L += y1[k] * w
    return L

x0, y0 = input_points()

a = float(input("请输入一个数值: "))
s = lagrange(a, x0, y0)
print("在输入点处的函数值为:", s)

t = pow(a,0.5)
#以计算x的开方为例 x0=100 x1=121 x2=144
#y0=10 y1=11 y2=12
#测试点x0=125 y0=11.17391
print("在测试点处的误差为:", s-t)

# 绘制插值结果图形
x_values = np.linspace(min(x0), max(x0), 1000)
y_values = [lagrange(x, x0, y0) for x in x_values]

plt.plot(x_values, y_values, label='Interpolation')
plt.scatter(x0, y0, color='red', label='Input Points')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Lagrange Interpolation')
plt.legend()
plt.grid(True)
plt.show()