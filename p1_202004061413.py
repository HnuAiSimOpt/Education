#引入matplotlib\numpy\random 库
import matplotlib.pyplot as plt
import numpy as np
import random

#定义ols函数
def ordinary_least_square( x0,list2):
    #计算下、y的平均值
    x_mean = sum(x0)/len(x0)
    y_mean = sum(y0)/len(y0)
    #计算平方和    
    x_square_sum = sum(x**2 for x in x0) 
    x_sum_square = (sum(x for x in x0))**2
   
    mul = [m*n for m,n in zip(x0,y0)]
    mul_sum = sum(mul)
     
    # 用公式求解线性方程的系数
    a = (len(x0)*mul_sum-sum(x0)*sum(y0))/(len(x0)*x_square_sum-x_sum_square)
    b = y_mean - a*x_mean

    # 拟合直线开始和结束的y值
    line_str_y = a*x0[0] + b
    line_end_y = a*x0[-1] + b

    fit_x = [x0[0], x0[-1]]
    fit_y = [line_str_y,line_end_y]
    #作图
    plt.scatter(x0,y0)
    plt.plot(fit_x,fit_y)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()
    print('y=',a,'x+',b)
    #生成有线性关系的数据
if __name__ == "__main__":
    x0 =np.arange(1,101)
    y0 =np.arange(1,101,1)
    
    for i in range(100):
    
        y0[i]=y0[i]+random.randrange(-10,20)
    #调用ols函数    
    ordinary_least_square( x0,y0)