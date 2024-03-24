#杨正昱 202104060905



def lagrange(x_known, y_known, x_new):
    n = len(x_known)
    y_new = 0
    for i in range(n):
        p = y_known[i]
        for j in range(n):
            if j != i:
                p *= (x_new - x_known[j]) / (x_known[i] - x_known[j])
        y_new += p
    return y_new


import random
x_known=sorted(random.sample(range(1,20),10))
y_known=random.sample(range(1,20),10)
x_new=eval(input("please input \nx_new="))
y_new=lagrange(x_known, y_known, x_new)
print("y_new=",y_new)




# 访问 https://www.jetbrains.com/help/pycharm/ 获取 PyCharm 帮助
