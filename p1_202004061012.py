import matplotlib.pyplot as plt
from pylab import mpl
"""一元线性拟合
@author  :   徐万青 202004061012 车辆2004班
"""
"""输入"""


"""
生成一组点集，设置初始函数为y=2x+3
"""
x = []
y = []
for i in range(20):
    x.append(i)
    a=2*i+3
    y.append(a)
print(x,y)

size = len(x)

"""拟合曲线参数计算"""
def nihe(data_x,data_y):
      i=0
      sum_x=0   
      sum_y=0
      sum_xy=0
      sum_sqare_x=0
      average_x=0
      average_y=0
      while i<size:
          sum_x+=data_x[i]                #x和计算
          sum_y+=data_y[i]                #y和计算
          sum_xy+=data_x[i]*data_y[i]     #xy和计算
          sum_sqare_x+=data_x[i]*data_x[i] #x方和计算
          i+=1
      average_x=sum_x/size                #x均值计算
      average_y=sum_y/size                #y均值计算
      a=(size*sum_xy-sum_x*sum_y)/(size*sum_sqare_x-sum_x*sum_x)#斜率a计算
      b=average_y-average_x*a             #截距b计算
      print("线性回归方程为y=",a,"*x+",b)  #输出拟合函数
      return [a,b]

"""拟合后曲线上相应的函数值的计算"""
def shiji(data_x,a,b):
    dy=[]
    for x in data_x:
        dy.append(a*x+b)               #求出拟合值并存入dy数组
    return dy

parameter = nihe(x,y)                     #调用拟合函数进行拟合

"""计算均方根误差"""
m=0
n=0
for m in range(size):
    n+=(shiji(x,parameter[0],parameter[1])[m]-y[m])**2#平方误差计算
    m+=m
print("均方根误差为",(n/size)**0.5)        #输出均方根误差
        
"""绘制函数"""
def draw(data_x,data_y_new,data_y_old):
    plt.plot(data_x,data_y_new,label="拟合曲线",color="black")
    plt.scatter(data_x,data_y_old,label="离散数据")
    mpl.rcParams['font.sans-serif'] = ['SimHei']
    mpl.rcParams['axes.unicode_minus'] = False
    plt.title("一元线性拟合数据")
    plt.legend(loc="upper left")
    plt.show()
 
 

draw_data = shiji(x,parameter[0],parameter[1])
draw(x,draw_data,y)
