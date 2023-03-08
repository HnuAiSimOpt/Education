# Test
For education
import ast
n=int(input("请输入数据的个数:"))
a=n-1#a为列表索引的下标
x=ast.literal_eval(input("请输入x的值，中间用逗号隔开:"))
y=ast.literal_eval(input("请输入y的值，中间用逗号隔开:"))
#键盘输入x，y的数据并转化成列表
def average(x):
    x_sum=0
    for i in range(0,a+1):
        x_sum +=x[i]
    x_average=x_sum/n
    return x_average
#函数average用来算x,y的平均数
 
x_average=average(x)
y_average=average(y)
 
product_xy = 0
for i in range(0,a+1):
    product_xy += x[i]*y[i]
    #计算x[i]*y[i]的和
 
product_xx=0
for i in range(0,a+1):
    product_xx += x[i]*x[i]
    #计算x[i]平方之和
 
b_up=product_xy - n*x_average*y_average#斜率分子的值
b_down=product_xx-n*(x_average**2)#斜率分母的值
b=b_up/b_down#斜率的值
 
a=y_average-b*x_average
print("截距的值为:%10f"%a)
print("斜率的值为:%10f"%b)
