#崔宗文 202104060926 车辆2102
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['font.family'] = ['sans-serif']    # 显示中文标签
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

def lagrange(x0, y0, x):
    n = x0.shape[0]
    s = 0
    for i in range(0, n):
        k = 1

        for j in range(0, n):
            if i != j:
                k = k*(x-x0[j])/(x0[i]-x0[j])
        s += k*y0[i]
    a = s
    return a

def fun(x):
    return 1/(1+x**2)

x0 = np.linspace(-5, 5, 10)   
y0 = fun(x0)   

x1 = np.linspace(-5, 5, 100)  # 为画图取得的100个点
y1 = lagrange(x0, y0, x1)     # 将100个点带入10插值多项式得出对应的函数值
f1 = lagrange(x0, y0, 3.5)    # 3.5处的插值
f2 = lagrange(x0, y0, 4.5)    # 4.5处的插值

plt.figure()
plt.title('lagrange插值多项式')
plt.plot(x1, fun(x1))   # 真实函数的曲线
plt.plot(x1, y1, 'r')   # 插值曲线
plt.plot(x0, y0, 'b.')  # 插值点
plt.plot(3.5, f1, 'g+')  # 两个特殊值点
plt.plot(4.5, f2, 'k+')

plt.legend(['真实函数值', '插值曲线', '实际取值点', '3.5处的插值函数值', '4.5处的插值函数值'])
plt.show()

print('3.5插值结果：', f1, '实际结果：', fun(3.5))
print('4.5插值结果：', f2, '实际结果：', fun(4.5))

