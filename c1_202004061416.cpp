#include <stdio.h>
#include <math.h>

#define N 10 // 样本数目

double x[N] = { 1,2,3,4,5,6,7,8,9,10 }; // 自变量数组
double y[N] = { 2.3,4.5,6.1,8.0,9.8,11.1,13.4,15.3,16.8,18.5 }; // 因变量数组

int main()
{
    double sumX = 0; // 自变量x的和
    double sumY = 0; // 因变量y的和
    double sumXY = 0; // 自变量x和因变量y乘积的和
    double sumX2 = 0; // 自变量x平方和

    for (int i = 0; i < N; i++)
    {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2 += pow(x[i], 2); // 计算平方函数，需要包含math库
    }

    double b = (N * sumXY - sumX * sumY) / (N * sumX2 - pow(sumX, 2)); // 求回归系数b
    double a = (sumY - b * sumX) / N; // 求截距a

    printf("y = %.2lfx + %.2lf", b, a); // 输出回归方程

    return 0;
}
