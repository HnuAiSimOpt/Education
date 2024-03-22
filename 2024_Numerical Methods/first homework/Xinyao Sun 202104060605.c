/******************************************************************************

Welcome to GDB Online.
  GDB online is an online compiler and debugger tool for C, C++, Python, PHP, Ruby, 
  C#, OCaml, VB, Perl, Swift, Prolog, Javascript, Pascal, COBOL, HTML, CSS, JS
  Code, Compile, Run and Debug online from anywhere in world.

*******************************************************************************/
#include <stdio.h>//车辆2102班 孙鑫耀
#include <math.h>
int main()
{
    //定义变量
    double a;
    double b;
    double c;
    double d;
    double h;
    //定义函数
    double f(double e)
    {
        double y = sqrt(e);//可在此处修改函数
        return y;
    }
    //输入插值点
    printf("请输入第一个插值点");
    scanf("%le",&a);
    printf("请输入第二个插值点");
    scanf("%le",&b);
    printf("请输入第三个插值点");
    scanf("%le",&c);
    printf("请输入要计算的点");
    scanf("%le",&d);
    double la(double x)
    {
        double g = ((x-b)*(x-c)*f(a))/((a-b)*(a-c))+((x-a)*(x-c)*f(b))/((b-a)*(b-c))+((x-a)*(x-b)*f(c))/((c-a)*(c-b));
        return g;
    }
    h=la(d);
    printf("结果为%f",h);
    return 0;
}

