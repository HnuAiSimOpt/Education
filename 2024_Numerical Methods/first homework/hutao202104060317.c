#include <stdio.h>
#include <stdlib.h>
#include <time.h>

float lagrange(int x, int n, int x1[], int x2[]);

int main() {
    int n, i, Xceshidian, Yceshidian;
    float Ylilunzhi, delta;

    printf("please input n: ");
    scanf("%d", &n);//输入想要的插值点组数

    int *X1 = (int *)malloc(n * sizeof(int));
    int *Y1 = (int *)malloc(n * sizeof(int));//给X1，Y1分配内存

    srand((unsigned int)time(NULL));//以当前时间为种子，以便后续随机数生成横、纵坐标的不一样
    for (i = 0; i < n; i++) {
        X1[i] = rand() % 100;
        printf("x%d\n", X1[i]);
    }//生成不同的插值点的横坐标

    for (i = 0; i < n; i++) {
        Y1[i] = rand() % 100;
        printf("y%d\n", Y1[i]);
    }//生成不同的插值点的纵坐标

    Xceshidian = rand() % 100;
    Yceshidian = rand() % 100;//生成一组测试点
    printf("%d\n%d\n", Xceshidian, Yceshidian);

    Ylilunzhi = lagrange(Xceshidian, n, X1, Y1);//调用lagrange函数
    delta = Yceshidian - Ylilunzhi;
    printf("Delta: %f\n", delta);//输出误差

    free(X1);
    free(Y1);//释放内存
    
    return 0;
}
//下面是lagerange函数
float lagrange(int x, int n, int x1[], int x2[]) {
    int i, j;
    float l1, s1, s2;

    s2 = 0;
    for (i = 0; i < n; i++) {
        s1 = 1;
        for (j = 0; j < n; j++) {
            if (i == j) {
                continue;
            } 
            else {
                if (x1[i]!=x1[j]){
                    l1 = (float)(x - x1[j]) / (x1[i] - x1[j]);
                    s1 = s1 * l1;
                }
                else{
                    printf("请重新输入n");
                    break;//当出现横坐标值相同时，终止程序，重新输入n，重新生成新的插值点
                }
            }
        }
        s2 = s2 + x2[i] * s1;
    }
    return s2;
}
