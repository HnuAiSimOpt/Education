#车辆2102班 陈奕瑞 202104060904
import numpy as np

#定义拉格朗日函数
def lagrange_interpolation(x_old, y_old, x_new):
    n = len(x_old)
    y_new = 0
    for i in range(n):
        p = y_old[i]
        for j in range(n):
            if j != i:
                p *= (x_new - x_old[j]) / (x_old[i] - x_old[j])
        y_new += p
    return y_new

# 已知的插值节点
x_old = np.array([0,1,2,3,4,5])
y_old= np.array([0,2.5,5,8,13.5,22])

#插值点
x_new = [1.2,2.5,3.6,4.3]

# 计算各个插值点的函数值
y_new = lagrange_interpolation(x_old, y_old, x_new)
for i, j in zip(x_new, y_new):
    print(f"插值点 x = {i} 的函数值 y = {j}")