#何子麒202104060408

def lagrange(x, y, num_points, x_test):
    # 所有的基函数值，每个元素代表一个基函数的值
    l = np.zeros(shape=(num_points, ))
    for k in range(num_points):
        l[k] = 1
        for k_ in range(num_points):
            if k != k_:
                l[k] = l[k]*(x_test-x[k_])/(x[k]-x[k_])
            else:
                pass


    L = 0
    for i in range(num_points):
        L += y[i]*l[i]
    return L