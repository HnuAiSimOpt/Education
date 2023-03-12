#include <stdio.h>

int main()
{
    float x[5] = {1, 3, 4, 5, 8};//输入X值
    float y[5] = {4, 2, 6, 4, 9}; //输入Y值
    float sum_x = 0, sum_y = 0;
    float sum_xx = 0, sum_xy = 0;
    float k, b;
    float e_2 = 0;
    int n = 5;  // 数据点的数目
    

    // 计算均值
       for (int i = 0; i < n; i++) {
        sum_x += x[i];
        sum_y += y[i];
    }
    float ave_x = sum_x / n;
    float ave_y = sum_y / n;


     // 计算X,y分别与平均值差和的乘积、x与平均值的差平方和
    for (int i = 0; i < n; i++) {
        sum_xy += (x[i] - ave_x) * (y[i] - ave_y);
        sum_xx += (x[i] - ave_x) * (x[i] - ave_x);
        
    }

    // 计算斜率b和截距a
    k = sum_xy / sum_xx;
    b = ave_y - k * ave_x;

    printf("斜率k为：%.2f\n", k);
    printf("截距b为：%.2f\n", b);
    printf("拟合直线方程为：y = %.2fx + %.2f\n",k ,b);
    
    //误差分析
    for (int i = 0; i < n; i++) {
        e_2 += (k * x[i] + b - y[i]) * (k * x[i] + b - y[i]);//误差平方和
        
    }
    printf("误差平方和e_2为：%.2f\n", e_2);
    return  0;
}