//徐雕锐 202104060420
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 10
double lagrange(double x[],double y[],int n,double t);
int main()
{
    double x[N]={0,1,2,3,4,5,6,7,8,9};
    double y[N]={0,1,4,9,16,25,36,49,64,81};
    double t;
    scanf("%lf",&t);
    printf("%lf\n",lagrange(x,y,N,t));
    return 0;
}
double lagrange(double x[],double y[],int n,double t)
{
    double s=0;
    for(int i=0;i<n;i++)
    {
        double p=1;
        for(int j=0;j<n;j++)
        {
            if(j!=i)
            {
                p*=(t-x[j])/(x[i]-x[j]);
            }
        }
        s+=y[i]*p;
    }
    return s;
}
