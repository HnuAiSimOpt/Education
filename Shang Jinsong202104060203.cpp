// ConsoleApplication2.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>

//尚锦松202104060203

// 求解拉格朗日插值
double lagrangeInterpolation(double x[], double y[], int n, double xi)
{
    double result = 0.0;

    for (int i = 0; i < n; i++)
    {
        double term = y[i];
        for (int j = 0; j < n; j++)
        {
            if (j != i)
            {
                term *= (xi - x[j]) / (x[i] - x[j]);
            }
        }
        result += term;
    }

    return result;
}

int main()
{
    // 给定的数据点
    double x[] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
    double y[] = { 1.0, 2.0, 11.0, 33.0, 65.0 };
    int n = sizeof(x) / sizeof(x[0]);

    // 要进行插值的点
    double xi = 3.0;

    // 求解插值
    double interpolatedValue = lagrangeInterpolation(x, y, n, xi);

    printf("在 x = %.2f 处的插值结果为: %.2f\n", xi, interpolatedValue);

    return 0;
}