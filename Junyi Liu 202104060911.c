#include <stdio.h>
//刘俊義202104060911
double lagrange_interpolation(double x[], double y[], int n, double xi) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        double term = y[i];
        for (int j = 0; j < n; j++) {
            if (i != j) {
                term *= (xi - x[j]) / (x[i] - x[j]);
            }
        }
        result += term;
    }
    return result;
}

int main() {
    int n;
    printf("输入数据点的数量: ");
    scanf("%d", &n);

    double x[n], y[n];
    printf("输入数据点 (x, y):\n");
    for (int i = 0; i < n; i++) {
        printf("x%d, y%d: ", i, i);
        scanf("%lf %lf", &x[i], &y[i]);
    }

    double xi;
    printf("输入要插值的点: ");
    scanf("%lf", &xi);

    double interpolated_value = lagrange_interpolation(x, y, n, xi);
    printf("插值位于 x=%.2f is y=%.2f\n", xi, interpolated_value);

    return 0;
}