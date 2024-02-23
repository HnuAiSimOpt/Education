#include <stdio.h>
#include <math.h>
// 定义一个表示点的结构体
struc
s
struct Point {    
    fl
float x, y;
};
// 计算两点之间的距离
float distance(struct Point p1, struct Point p2) {
    
    ret
return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

// 计算三角形的周长
flo
float calculatePerimeter(struct Point p1, struct Point p2, struct Point p3) {
    
    flo
 
float side1 = distance(p1, p2);
float side2 = distance(p2, p3);      
float side3 = distance(p3, p1);
    ret
return side1 + side2 + side3;
}
}
// 计算三角形的面积
f
float calculateArea(struct Point p1, struct Point p2, struct Point p3) {    
    floa
float side1 = distance(p1, p2);      
float side2 = distance(p2, p3);   
    floa
float side3 = distance(p3, p1);
// 使用海伦公式计算面积
    float s = (side1 + side2 + side3) / 2;  
return sqrt(s * (s - side1) * (s - side2) * (s - side3));
}
}
i
int main() {   
    str
struct Point p1, p2, p3;
    // 输入三个点的坐标   
printf("Enter coordinates for point 1 (x y): ");   
scanf("%f %f", &p1.x, &p1.y);
    printf("Enter coordinates for point 2 (x y): ");
    scanf("%f %f", &p2.x, &p2.y);
    prin
printf("Enter coordinates for point 3 (x y): ");
    scanf("%f %f", &p3.x, &p3.y);
    // 计算并输出周长和面积
float perimeter = calculatePerimeter(p1, p2, p3);
    float area = calculateArea(p1, p2, p3);
    pr
printf("Perimeter of the triangle: %.2f\n", perimeter);  
    pri
printf("Area of the triangle: %.2f\n", area);
    retu
return 0;
}
