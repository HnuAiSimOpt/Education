/*许少辉 
车辆2003班
202004061301 
2023.03.06
通过公式直接求回归方程，在该代码下面有通过梯度下降来求损失函数的最小值的另一端代码，但代码拟合误差太大且未找到原因，
我用matlab用同样的方法写了一遍，奇怪的是用matlab用同样的算法写的就能很好的拟合，同时提交了matlab的文件望老师能够指出代码错误 


*/ 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 5/* 定义数据点数 */
double X[N] = {1, 2, 3, 4, 5};/* 定义X和Y数组，存储数据点 */
double Y[N] = {1.2, 2.5, 3.6, 4.8, 6.1};
int main()
{
    double sum_x = 0, sum_y = 0; /* X和Y的和 */
    double sum_xy = 0, sum_x2 = 0; /* XY和X^2的和 */
    double a, b; /* 拟合直线的系数 */
    for (int i = 0; i < N; i++)/* 计算X、Y、XY、X^2的和 */
    {
        sum_x += X[i];
        sum_y += Y[i];
        sum_xy += X[i] * Y[i];
        sum_x2 += X[i] * X[i];
    }
    b = (N * sum_xy - sum_x * sum_y) / (N * sum_x2 - pow(sum_x, 2));/* 计算拟合直线的系数 */
    a = (sum_y - b * sum_x) / N;
	printf("拟合直线方程为: y = %.2fx + %.2f", b, a);
    return 0;
}
/*
#include <iostream>
#include <stdio.h>
#include<math.h>
#define v 0.000001   //设置下山步长 
int main()
{
    using namespace std;
    int i,j,k,m;
    int col;
    int number;
    double sum_x = 0;
    double sum_y = 0;
    double  x[100];         //数据样本
    double  y[100];
    double theta1;       //一次项的系数
    double theta0;       //常系数 
    double e;            //误差系数
    double old_theta1;
    double old_theta0;
    printf("please input sample's number :"); //输入样本的个数
    cin>>number;
    printf("sample's number is %d\n\n",number);
    printf("please input sample'x :");       //样本的输入x
    for(k=0;k<number;k++)
    {
        cin>>x[k];
    }
    printf("\n");
    printf("please input sample'y :");      //输入的样本y
    for(m=0;m<number;m++)
    {
        cin>>y[m];
    }
    printf("\n");
    col = number;
    printf("\n\n");
    printf("the training sample is : %d \n",col);
    for(i=0;i<col-1;i++)
    {
        sum_x = sum_x + x[i];
        sum_y = sum_y + y[i];
    }
    theta1 = (double)sum_y/sum_x;//设置下山的初始点
    theta0 = y[3]-theta1*x[3];
    while(1)                     //梯度下降，直到找到最优解退出循环
    {
       double temp1 = 0;
       double temp0 = 0;
       for(j=0;j<col-1;j++)//计算损失函数分别对常系数和一次项系数的导数
       {
          temp1 = temp1 + ((theta0 + theta1*x[j])-y[j])*x[j];
          temp0 = temp0 + ((theta0 + theta1*x[j])-y[j])*1;
       }
	   temp1 = temp1 / col;
       temp0 = temp0 / col;
       old_theta1 = theta1;//将前一个常系数和一次项系数存储以后续比较
       old_theta0 = theta0;

       theta1 = theta1 - v*temp1;//更新每个样本的常系数和一次项系数
       theta0 = theta0 - v*temp0;
       temp1 = 0;
       temp0 = 0;
       e = (pow((old_theta1-theta1),2) + pow((old_theta0 - theta0),2));
        if(e<0.000003);//设置下降速度的一个阈值
        {
            printf("the objective function is : \n");
            printf("f(x)= %f + %f*x",theta0,theta1);
            break;//跳出循环
        }
    }
}
*/
