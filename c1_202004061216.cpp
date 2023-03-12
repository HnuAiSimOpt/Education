/* 
1.目的：使用最小二乘法对数据点进行线性拟合
2.作者：高旭根
  日期：2023年3月8日
  学号：202004061216
3.输入：给定输入
  输出：拟合函数表达式
        均方根误差
*/ 


#include <stdio.h>
#include "time.h"
#include <stdlib.h>
#include <math.h>


int main()  {
	float x[4] = {2.1, 4.0, 7.5, 9.4} ;
	float y[4] = {7.3, 11.3, 17.2, 21.7};//参照y=2x+3给出一组数据
	float y_test[4];
	float average_x, average_y, sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
	float k, b;//k为拟合函数的截距，b为拟合函数的斜率；
	float error, error_sum = 0.0;
	int i;
	
	
    for(i = 0; i < 4; i++) {         
	    sum_x = sum_x + x[i];//求x的和  
		sum_y = sum_y + y[i];//求y的和 
		sum_xy = sum_xy + x[i] * y[i];//求x*y的和 
		sum_x2 = sum_x2 + x[i] * x[i];//求x^2的和
	}
    average_x = sum_x / 4;//求x的均值 
    average_y = sum_y / 4;//求y的均值 
    k = (sum_xy - 4 * average_x * average_y) / (sum_x2 - 4 * average_x * average_x);//最小二乘法求拟合函数斜率 
    b = average_y - k * average_x;//所求为拟合函数截距
    printf("拟合函数为：y = %fx + %f\n", k, b);
    
    //计算均方根方差 
	for(i = 0; i < 4; i++) {	
		y_test[i] = k * x[i] + b;
		error = fabs(y_test[i] - y[i]);//求得拟合函数与实际函数函数值的误差 
		error_sum = error_sum + error * error;//求出该误差的平方和 
	    }
	printf("均方根误差为：%f", sqrt(error_sum / 4));//求出均方根误差，并输出 
    return 0;
}    
