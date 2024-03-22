/*

罗珺  202104060809

*/


#include <stdio.h>

// 定义拉格朗日插值函数
float lagrange(float x, int n, float X[], float Y[]) {
    float result = 0.0;  // 初始化插值结果

    for (int i = 0; i < n; i++) {
        float lk = Y[i];  
        
        // 计算每个基函数的乘积 y*lk(x)
        for (int j = 0; j < n; j++) {
            if (j != i) {
                lk *= (x - X[j]) / (X[i] - X[j]);
            }
        }

        // 基函数的乘积相加得到插值结果
        result += lk;
    }

    return result;
}
int main() {
    int n;
    printf("输入数据点个数: ");
    scanf("%d", &n);

    float X[n], Y[n];

    // 输入数据点
    for (int i = 0; i < n; i++) {
        printf("输入第%d个数据点X值: ", i+1);
        scanf("%f", &X[i]);

        printf("输入第%d个数据点Y值: ",i+1);
        scanf("%f", &Y[i]);
    }

    float x;
    printf("需要计算插值点的X值: ");
    scanf("%f", &x);

    // 进行拉格朗日插值
    float result = lagrange(x, n, X, Y);
    printf("插值点：(%.2f , %.2f)\n", x, result);

    return 0;
}
