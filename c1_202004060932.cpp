/******************************************************************************
author:段蕴 
class：车辆2004班
date：2023/3/8
Student ID：202004060932
theme：OLS
purpose：给定范围生成2组数x和y，每组随机生成20个数并从小到大排列，依次对应，应用OLS最小二乘法求解回归方程y=b0x+b1
(*^_^*)：图像暂时还不会用C++画，还望助教指点，谢谢。
*******************************************************************************/
#include <stdio.h>
#include <string.h>  
#include <math.h> 
#include <time.h> 
#include<cstdlib>

#define NR(x) sizeof(x)/sizeof(x[0])  
float zxec(const double *x, const double *y, int n, double *b1, double *b0)//定义最小二乘函数各个参数值   
{  
    int i;  
    double  sx,sy,sx2,sxy,error;  
    sx = 0.0;  
    sy = 0.0;  
    sx2 = 0.0;  
    sxy = 0.0;
    error = 0.0;
    for (i = 0; i < n; i++) {  
      sx = sx + x[i];//计算横坐标X的和  
      sy = sy + y[i];//计算纵坐标y的和   
      sx2 = sx2 + pow(x[i], 2.0);//计算x的平方和  
      sxy = sxy + (x[i] * y[i]);//计算x*y的和  
    }  
    //根据公式求解b1和b0的值   
    *b1 = (sxy - ((sx * sy)/(double)n)) / (sx2-(pow(sx,2.0)/(double)n));  
    *b0 = (sy - ((*b1) * sx)) / (double)n;
    for (int i = 0; i < n; i++){
        error = error + (y[i] - (*b1 * x[i] + *b0)) * (y[i] - (*b1* x[i] + *b0)); //计算拟合误差
    }
    return error;
}

int main()
{
    //1~20随机生成20个x值
    int i;
    double x[20]; 
    printf("20 x numbers：\n");
    srand(time(0)); //保证每个数字大概率不相同
    for (i=0; i<20; i++)
    {
        x[i]=rand()%20+1;
    }
    //x排序
    int iTx;//定义变量，表示最小的数组元素 
	int iPx;//定义变量，表示元素位置 
	int j;
    //使用选择法对数组元素从小到大排序
	for(i=0;i<19;i++)//外层循环下标为0-18，表示前19个数字 
	{
		iTx=x[i];//假设当前数字为最小值 
		iPx=i;//记录最小元素位置 
		for(j=i+1;j<20;j++)//设置内层循环下标为i+1-19，表示剩下的未排序数组元素 
		{
			if(x[j]<iTx)//如果后续元素中有比前面设定的最小值还小 
			{
				iTx=x[j];//重新设定最小值 
				iPx=j;//修正最小元素位置 
			}
		}
		x[iPx]=x[i];//将最小的数组元素和当前排序次数对应的数组元素互换 
		x[i]=iTx;
	 } 
	 for(i=0;i<20;i++)//输出x数组 
	 {
	 	printf("%10f",x[i]);
	 }
	 printf("\n");
	 
    //同理1~60随机生成20个y值
    int k;
    double y[50];  
    printf("\n20 y numbers :\n");
    srand(time(0));  
    for (k=0; k<20; k++)
    {
        y[k]=rand()%60+1; 
    }
    //y排序
    int iTy; 
	int iPy;
	int l;
	for(i=0;i<19;i++)
	{
		iTy=y[i]; 
		iPy=i; 
		for(l=i+1;l<20;l++) 
		{
			if(y[l]<iTy) 
			{
				iTy=y[l]; 
				iPy=l; 
			}
		}
		y[iPy]=y[i]; 
		y[i]=iTy;
	 } 
	 for(i=0;i<20;i++)//输出y数组 
	 {
	 	printf("%12f",y[i]);
	 }
	 printf("\n");
     //线性拟合
	 double b1=0,b0=0;
	 printf("误差大约为：%f\n", zxec(x, y, 20, &b1, &b0));
	 printf("上述两组数据相应的回归方程为:y=%fx+%f\n", b1,b0);//输出拟合方程
}


