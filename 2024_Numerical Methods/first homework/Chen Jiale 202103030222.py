#车辆2102班 陈嘉乐 202103030222
def lagrange(a, b, x):
    result = 0.0
    n = len(a)
    for i in range(n):
        term = b[i]
        for j in range(n):
            if j != i:
                term *= (x - a[j]) / (a[i] - a[j])
        result += term
    return result
#所有数据点
a = input("x: ").split()
b = input("y: ").split()
a = list(map(float, a))
b = list(map(float, b))
# 根据索引列表提取需要插值的数据点
index_list = input("插值点对应的索引 ").split()
index_list = list(map(int, index_list))
x_data = [a[i] for i in index_list]
y_data = [b[i] for i in index_list]
#插值
def poly(x):
    return lagrange(x_data, y_data, x)

# 评估插值多项式1
x = float(input("插值曲线上任意一点x: "))
result = poly(x)
print(f"x的纵坐标: {result}")
