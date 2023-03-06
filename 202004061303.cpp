#include <stdio.h>
//定义最小二乘法函数
float LSquareFitting(float* x, float* y, int n, float* a, float* b)
{
	float sumx = 0.0;
	float sumy = 0.0;
	float sumx2 = 0.0;
	float sumxy = 0.0;
	float errorSquare = 0;
	for (int i = 0; i < n; i++)
	{
		sumx += x[i];
		sumy += y[i];
		sumx2 += x[i] * x[i];
		sumxy += x[i] * y[i];
	}
	*a = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx * sumx);//最小二乘法一次项系数a
	*b = ((sumx2 * sumy) - (sumx * sumxy)) / (n * sumx2 - sumx * sumx);//最小二乘法常数项系数b
	for (int i = 0; i < n; i++)
	{
		errorSquare = errorSquare + (y[i] - (*a * x[i] + *b)) * (y[i] - (*a * x[i] + *b)); //最小二乘法拟合误差
	}
	return errorSquare;
}


int main()
{
	float aa, bb;
	float x[30] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30};
	float y[30] = { 4,7,12,19,28,34,40,48,57,67,78,89,96,103,112,119,127,134,143,150,157,170,179,185,192,200,210,218,226,235};
	//调用函数并输出拟合误差和拟合函数
	printf("%f\n", LSquareFitting(x, y, 30, &aa, &bb));
	printf("y=%fx%f\n", aa,bb);
	return 0;
}