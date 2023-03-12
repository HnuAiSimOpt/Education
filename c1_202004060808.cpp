/* 
1.目的：使用最小二乘法对数据点进行线性拟合
2.作者：罗兆锋
  日期：2023年3月8日
  学号：202004060808
3.输入：随机数输入
  输出：拟合函数表达式
        均方根误差
*/ 



#include <stdio.h>
#include "time.h"
#include <stdlib.h>
#include <math.h>

int main() //实际函数为 y = x + 1,通过在实际点上随机施加误差，得到拟合所需数据点，后对这些点做一次线性拟合 
{
	int size = 100, size_test = 1000;//拟合所用数字数量为100，检验函数拟合情况所用数字数量为1000 
    float x[100], y[100];//第一个数组存储0-100之间的随机数，第二个数组存储第一个数组每个数在某个二次函数的实际函数值 
	float x_test[1000], y_test[1000]; //第一个数组存储0-100之间的随机数，第二个数组存储第一个数组每个数在某个二次函数的实际函数值 
    float a, b;//a为拟合函数的截距，b为拟合函数的斜率； 
	float average_x, average_y, sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
	float error, error_test, error_test_sum = 0.0;
    int i, j;

    srand(time(0));
    for(i=0; i<size; i++) {        
	    x[i] = 0+1.0*(rand()%RAND_MAX)/RAND_MAX *(100-0);//随机在[0,100]产生一个x值 
	    error = (2.0 * rand() / RAND_MAX - 1.0) * 0.1;//随机产生10%的误差 
	    sum_x = sum_x + x[i];//求x的和 
	    y[i] = (x[i] + 1) * (1.0 + error);//令真实函数值上下有10%的跳动，即产生一个小的误差 
		sum_y = sum_y + y[i];//求y的和 
		sum_xy = sum_xy + x[i] * y[i];//求x*y的和 
		sum_x2 = sum_x2 + x[i] * x[i];//求x^2的和 
	}
    average_x = sum_x / size;//求x的均值 
    average_y = sum_y / size; //求y的均值 
    b = (sum_xy - size * average_x * average_y) / (sum_x2 - size * average_x * average_x);//最小二乘法求拟合函数斜率 
    a = average_y - b * average_x;//所求为拟合函数截距 
    if(a >= 0) printf("拟合函数为：%fx+%f\n", b, a);
    if(a < 0) printf("拟合函数为：%fx%f", b, a);//输出拟合函数表达式 
    for(i=0; i<size_test; i++) {
        x_test[i] = 0+1.0*(rand()%RAND_MAX)/RAND_MAX *(100-0);
        for(j = 0; j < size; j++) {
        	if(x_test[i] == x[j]) {
        		i--;
        		break;
			    }
            }//令第二次所取x值与第一次不同 
		y_test[i] = b * x_test[i] + a;
		error_test = fabs(y_test[i] - (x_test[i] + 1));//求得拟合函数与实际函数函数值的误差的平方和 
		error_test_sum = error_test_sum + error_test * error_test;
	    }
	printf("均方根误差为：%f", sqrt(error_test_sum / size_test));//求出均方根误差，并输出 
    return 0;
}
