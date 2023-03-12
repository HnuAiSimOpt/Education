#include <stdio.h>

float zuoye(float* x,float* y,int n, float* a,float* b)
{
	float sumx=0.0;
	float sumy=0.0;
	float sumx2=0.0;
	float sumxy=0.0;
	
	for(int i=0;i<n;i++)   /*进行累加*/
	{
		sumx += x[i];
		sumy += y[i];
		sumx2 += x[i]*x[i];   
		sumxy += x[i]*y[i];     
	}
	*a = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx * sumx);   
	*b = ((sumx2 * sumy) - (sumx * sumxy))/(n * sumx2 - sumx * sumx);    /*由最小二乘法推得的公式*/
	
	
	
}
		 
int main()
{
   float a, b;    
	
	float x[] = {1,2,3,4,5};     /*此处输入横坐标*/
	float y[] = {2,5,8,9,11};  /*此处输入纵坐标*/
    
	
	zuoye(x,y,5, &a, &b);  /*运行子程序*/
	printf("y=%fx+(%f)\n",a,b); /*输出函数表达式*/
  
   
	
}