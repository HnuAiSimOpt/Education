/*
许少辉 202004061301 车辆2003班
最小二乘法来计算线性拟合 




*/
#include <iostream>
#include <stdio.h>
#include<math.h>
#define v 0.00001   /*设置下山步长 */
int main()
{
    using namespace std;
    int i,j,k,m;
    int col;
    int number;
    int sum_x = 0;
    int sum_y = 0;
    int  x[100];         /*数据样本*/
    int  y[100];
    double theta1;       /* 一次项的系数*/
    double theta0;       /* 常系数 */
    double e;            /*误差系数*/
    double old_theta1;
    double old_theta0;
    printf("please input sample's number :"); /*输入样本的个数*/
    cin>>number;
    printf("sample's number is %d\n\n",number);
    printf("please input sample'x :");       /*样本的输入x*/
    for(k=0;k<number;k++)
    {
        cin>>x[k];
    }
    printf("\n");
    printf("please input sample'y :");      /*输入的样本y*/
    for(m=0;m<number;m++)
    {
        cin>>y[m];
    }
    printf("\n");
    k=0;
    m=0;
    printf("x[]= ");                /*输出样本点*/ 
    for(k=0;k<number;k++)
    {
        printf("%3d",x[k]);
    }
    printf("\n\n");
    printf("y[] = ");
    for(m=0;m<number;m++)
    {
       printf("%3d",y[m]);
    }
    col = number;
    printf("\n\n");
    printf("the training sample is : %d \n",col);
    for(i=0;i<col-1;i++)
    {
        sum_x = sum_x + x[i];
        sum_y = sum_y + y[i];
    }
    theta1 = (double)sum_y/sum_x;/*设置下山的初始点*/
    theta0 = y[0]-theta1*x[0];
    while(1)                     /*牛顿开始下山，直到找到最优解退出循环*/
    {
       double temp1 = 0;
       double temp0 = 0;
       for(j=0;j<col-1;j++)/*计算损失函数分别对常系数和一次项系数的导数*/
       {
          temp1 = temp1 + (y[j]-(theta0 + theta1*x[j]))*x[j];
          temp0 = temp0 + (y[j]-(theta0 + theta1*x[j]))*1;
       }
	   temp1 = temp1 / col;
       temp0 = temp0 / col;
       old_theta1 = theta1;/*将前一个常系数和一次项系数存储以后续比较*/
       old_theta0 = theta0;

       theta1 = theta1 - v*temp1;/*更新每个样本的常系数和一次项系数*/
       theta0 = theta0 - v*temp0;
       e = (pow((old_theta1-theta1),2) + pow((old_theta0 - theta0),2));
        if(e<0.000003);/*设置下降速度的一个阈值*/ 
        {
            printf("the objective function is : \n");
            printf("f(x)= %f + %f*x",theta0,theta1);
            break;/*跳出循环*/
        }
    }
}
