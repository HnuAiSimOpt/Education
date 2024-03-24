//汪淼 202104060520
#include <stdio.h>

double lagrange(int n, double x[], double y[], double test_x) {
    double test_y = 0;
    
    for(int i = 0; i < n; i++) {
        double sum_y = y[i];
        for(int j = 0; j < n; j++) {
            if (i != j) {
                sum_y *= (test_x - x[j]) / (x[i] - x[j]);
            }
        }
        test_y += sum_y;
    }
    return test_y;
}

int main() {
    int n = 10;
    double x[] = {1, 2, 4, 6, 8, 10, 12, 14, 16, 18};
    double y[] = {3, 6, 7, 8, 10, 12, 15, 17, 19, 20};
    double test_x, test_y;

    printf("Please input test_x: ");
    scanf("%lf", &test_x);

    test_y = lagrange(n, x, y, test_x);
    printf("test_y = %lf\n", test_y);

    return 0;
}
