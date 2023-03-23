//202004061311
//车辆2003 赵乙鑫
//2023.3.8
/*通过数组距离平方求偏导为0时有最佳拟合方案，得到a,b的求解公式，作为最小二乘的运算结果。包含输入数组的部分，可实现任意数量的数组进行拟合*/
/*其中
*a =(n*sumxy-sumx*sumy)/(n*sumx2-sumx*sumx);
*b=((sumx2*sumy)-(sumx*sumxy))/(n*sumx2-sumx*sumx);*/
#include<stdio.h>
//定义OLS函数
int LS_fun(float* x,float* y,int n,float*a,float*b,float*Error_Range)
{
	float sumx=0.0;//定义x总和
	float sumy=0.0;//定义y总和
	float sumxy=0.0;//定义x*y的总和
	float sumx2=0.0;//定义x平方的总和
	*Error_Range=0.0;//定义误差
	for(int i=0;i<n;i++)
	{
		sumx+=x[i];//求解x总和
		sumy+=y[i];//求解y总和
		sumx2+=x[i]*x[i];//求解x平方的总和
		sumxy+=x[i]*y[i];//求解x*y的总和
	}	
	//输入a,b的求解公式
		*a =(n*sumxy-sumx*sumy)/(n*sumx2-sumx*sumx);
		*b=((sumx2*sumy)-(sumx*sumxy))/(n*sumx2-sumx*sumx);
		for(int i=0;i<n;i++)
		{
			*Error_Range +=(y[i]-(*a*x[i]+*b))*(y[i]-(*a*x[i]+*b));//计算误差
		}
		
		return 0;
}
int main()
{
	float aa,bb,Error_Range;//定义输出的值
	int length;//定义数组长度
	scanf("%d",&length);//手动输入数组的长度
	printf("输入数组的长度为%d\n",length);
	float x[length];//定义x数组
	float y[length];//定义y数组
	int i;
	//设置循环实现依次输入x y数组内的值
	for(i-0;i<length;i++)
	{
	printf("input x[%d] ",i);
	scanf("%f",x+i);
	printf("input y[%d] ",i);
	scanf("%f",y+i);
	}
	LS_fun(x,y,length,&aa,&bb,&Error_Range);
	printf("%f\n",aa);//输出a的值
	printf("%f\n",bb);//输出b的值
	printf("%f\n",Error_Range);//输出误差的值
}
