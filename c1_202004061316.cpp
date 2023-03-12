#include <stdio.h>
#include <stdlib.h>

void linear_regression(double *x, double *y, int n, double *m, double *b);

int main()
{
    // 定义输入数据
    double x[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double y[5] = {0.5, 2.5, 2.0, 4.0, 3.5};
    
    // 定义输出参数
    double m, b;
    
    // 调用线性回归函数
    linear_regression(x, y, 5, &m, &b);
    
    // 输出结果
    printf("slope: %f", m);
    printf("intercept: %f", b);
    
    return 0;
}
void linear_regression(double *x, double *y, int n, double *m, double *b)
{
    // 计算 X 和 Y 的平均值
    double sumX = 0;
    double sumY = 0;
    
    for (int i = 0; i < n; i++) {
        sumX += x[i];
        sumY += y[i];
    }
    
    double avgX = sumX / n;
    double avgY = sumY / n;
    
    // 计算斜率和截距
    double numerator = 0;
    double denominator = 0;
    
    for (int i = 0; i < n; i++) {
        numerator += (x[i] - avgX) * (y[i] - avgY);
        denominator += (x[i] - avgX) * (x[i] - avgX);
    }
    
    *m = numerator / denominator;
    *b = avgY - (*m) * avgX;
}