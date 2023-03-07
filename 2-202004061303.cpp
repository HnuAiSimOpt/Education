//向人强 车辆2003
//解决的问题：OLS拟合。
//分为两节，第一节是定义函数，第二节是调用并输出拟合方程
#include <stdio.h>
#define DEFAULT_EPS 1e-15
#include <cmath>
//定义最小二乘函数并初始化
float LeastSquareLinearFit(double* x, double* y, const int num, double* a, double* b)
{
    int i = 0;
    double denominator = 0.0;
    double sum_xsquared = 0.0;
    double sum_y = 0.0;
    double sum_x = 0.0;
    double sum_xy = 0.0;
    double errorSquare = 0;
    //构建循环体
    for (i = 0; i < num; ++i)
    {
        sum_xsquared += x[i] * x[i];
        sum_y += y[i];
        sum_x += x[i];
        sum_xy += x[i] * y[i];//用循环叠加求解系数用到的中间量
    }


    denominator = (num * sum_xsquared - sum_x * sum_x);//计算系数a时的分母
    if (fabs(denominator) <=( DEFAULT_EPS))//判定分母不为0
    {
        return 1;
    }
    *a = (num * sum_xy - sum_x * sum_y) / denominator;//计算一次项系数a
    *b = (sum_xsquared * sum_y - sum_x * sum_xy) / denominator;//计算常数项b
    for (int i = 0; i < num; i++)
    {
        errorSquare = errorSquare + (y[i] - (*a * x[i] + *b)) * (y[i] - (*a * x[i] + *b)); //最小二乘法拟合误差
    }
    return errorSquare;
}

    
int main()
{
   double aa=0, bb=0;//初始化aa,bb

   //输入数组数据
    double x[30] = { 2,4,6,8,10,12,14,15,17,18,20,22,24,26,29};
    double y[30] = { 4,9,13,15,21,24,28,31,34,35,39,43,47,51,57};
    printf("%f\n", LeastSquareLinearFit(x, y, 30, &aa, &bb));
    //调用函数并输出拟合误差和拟合函数
    printf("y=%fx+%f\n", aa, bb);//输出拟合方程
    return 0;
}

