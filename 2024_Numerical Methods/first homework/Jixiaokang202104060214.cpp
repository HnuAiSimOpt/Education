//纪小康202104060214
#include <iostream>
using namespace std;

const int n=10;
//插值节点数量(预设)

float lagrange(float arrx[n],float arry[n],int n,float x)
//x[n],y[n]是插值节点坐标,n是插值节点数量与循环次数，x是要求取的插值节点的x值
{
	float lag[n];
	//插值基函数
	float y=0;
	//插值节点y值
	float up,down;
	//基函数上下两项

	for(int i=0;i<n;i++)
	{
		up=1;
		down=1;
		for(int j=0;j<n;j++)
		{
			if(j==i)
			{
				continue;
			}
			up*=(x-arrx[j]);
			down*=(arrx[i]-arrx[j]);
		}
		lag[i]=up/down;
	}
	for(int k=0;k<n;k++)
	{
		y+=arry[k]*lag[k];
	}
	return y;
}

int main()
{
	float arrx[n],arry[n];
	int num;
		cout<<"输入插值节点个数:";
	cin>>num;

	for(int i=0;i<num;i++)
	{
		cout<<"第"<<i+1<<"个节点的x值:";
		cin>>arrx[i];
		cout<<"第"<<i+1<<"个节点的y值:";
		cin>>arry[i];
	}

	float x;
	cout<<"要求取的插值节点的x值:";
	cin>>x;
		float result =lagrange(arrx,arry,num,x);
	cout<<"插值结果:"<<result;

	return 0;
}