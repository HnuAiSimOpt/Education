#Least square method
#车辆2004 陈靖彦 202004061107
#date:2023.3.3
#ordinary least squares
#自设置数据Xi=10,20,30,40,50,60,70,80,90,100
#自设置数据Yi=62,68,75,81,89,95,102,108,115,122
x=[10,20,30,40,50,60,70,80,90,100];
y=[62,68,75,81,89,95,102,108,115,122];
#求出拟合曲线的斜率与截距 y=bx+a
def OLS(data_x,data_y):
    b=0;#拟合の斜率
    a=0;#拟合の截距
    length=len(data_x);
    sum_Square_x=0
    sum_x=0;
    sum_y=0;
    sum_xy=0;
    average_x=0;
    average_y=0;
    Square_average_x=0;
    sum_Square_x=0;
    i=0;
    while i<length:
        sum_xy+=data_x[i]*data_y[i];
        sum_x +=data_x[i];
        sum_y +=data_y[i];
        sum_Square_x+=data_x[i]*data_x[i];
        i+=1;
    average_x=sum_x/length;
    average_y=sum_y/length;
    Square_average_x=average_x*average_x;
    b=(sum_xy-length*average_x*average_y)/(sum_Square_x-length*Square_average_x);
    a=average_y-b*average_x;
    return[b,a];

#计算拟合后对应函数值
[Fitting_b,Fitting_a]=OLS(x,y);
bx=[item*Fitting_b for item in x];
Fitting_y=[];
for i in bx:
    Fitting_y.append(i+Fitting_a)

print('y=','%.2f'%Fitting_b,'*x+','%.2f'%Fitting_a)

#计算残差平方和以衡量拟合数据质量

Subtraction=[c-d for c,d in zip(y,Fitting_y)];#实际值与拟合值作差

i=0;
while i<len(y):
    RSS=0;
    RSS+=Subtraction[i]*Subtraction[i];
    i+=1;
print("残差平方和RSS=",end="");
print('%.8f'%RSS);

#绘制图像
import matplotlib.pyplot as plt
from pylab import mpl
def draw(data_x,data_y_new,data_y_old):
    plt.plot(data_x,data_y_new,label="拟合曲线",color="black");
    plt.scatter(data_x,data_y_old,label="离散数据");
    mpl.rcParams['font.sans-serif'] = ['SimHei'];
    mpl.rcParams['axes.unicode_minus'] = False;
    plt.title("一元线性拟合数据");
    plt.legend(loc="upper left");
    plt.show();
 

parameter = OLS(x,y);
draw_data = Fitting_y;
draw(x,draw_data,y); 

    