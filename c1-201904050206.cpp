#包括<stdio.h>

int main()
{
    浮点数 x[5] = {1， 3， 4， 5， 8};//输入X值
    浮点数 y[5] = {4， 2， 6， 4， 9}; //输入Y值
    浮点数 sum_x = 0，sum_y = 0;
    浮点数 sum_xx = 0，sum_xy = 0;
    浮点 k， b;
    浮点数e_2 = 0;
    整数 n = 5;数据点的数目
    

    计算均值
       for （int i = 0; i < n; i++） {
sum_x += x[i];
sum_y += y[i];
    }
    浮点ave_x = sum_x / n;
    浮点ave_y = sum_y / n;


     计算X，y分别与平均值差和的乘积、x与平均值的差平方和
    for （int i = 0; i < n; i++） {
sum_xy += （x[i] - ave_x） * （y[i] - ave_y）;
sum_xx += （x[i] - ave_x） * （x[i] - ave_x）;
        
    }

    计算斜率b和截距a
k = sum_xy / sum_xx;
b = ave_y - k * ave_x;

    printf（“斜率k为：%.2f\n”， k）;
    printf（“截距B为：%.2f\n”， b）;
    printf（“拟合直线方程为：y = %.2fx + %.2f\n”，k ，b）;
    
    //误差分析
    for （int i = 0; i < n; i++） {
e_2 += （k * x[i] + b - y[i]） * （k * x[i] + b - y[i]）;//误差平方和
        
    }
    printf（“误差平方和e_2为：%.2f\n”， e_2）;
    返回 0;
}
