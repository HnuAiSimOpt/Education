#车辆2101黄自强 202104060609

import numpy as np  
def Lagrange(data_x, data_y, x):
    y = 0
    for i in range(len(data_y)):
        term = data_y[i]
        for j in range(len(data_x)):
            if j != i:
                term = term*(x-data_x[j])/(data_x[i]-data_x[j])
        y=y+term
    return y
data_x= np.array(input("请输入数据点的x坐标，用空格分隔：").split(), dtype=float)  
data_y= np.array(input("请输入数据点的y坐标，用空格分隔：").split(), dtype=float)  
x = float(input("请输入需要插值的x坐标："))  
y = Lagrange(data_x, data_y, x)  
print("插值结果: y = {y}".formate(y=y))

