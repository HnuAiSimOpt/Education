# 2维线性三角形单元python语言求解，使用numpy库，输入示例见135行
import numpy as np


def element_triangle(E, miu, t, node_ele, ele_node, solve):
    # 用行列式计算常量
    arr = np.concatenate((np.ones([3, 1]), node_ele[ele_node, ...]), axis=1)
    # print(arr)
    b1 = -np.linalg.det(arr[[1, 2]][..., [0, 2]])
    # print(b1)
    b2 = -np.linalg.det(arr[[2, 0]][..., [0, 2]])
    b3 = -np.linalg.det(arr[[0, 1]][..., [0, 2]])
    c1 = np.linalg.det(arr[[1, 2]][..., [0, 1]])
    c2 = np.linalg.det(arr[[2, 0]][..., [0, 1]])
    c3 = np.linalg.det(arr[[0, 1]][..., [0, 1]])
    A = np.abs(np.linalg.det(arr) / 2)
    if A == 0:
        raise Exception("三角形单元面积为零，请检查单元节点序号{}是否有误".format(ele_node))
    # print(A)
    B = np.array([[b1, 0, b2, 0, b3, 0],
                  [0, c1, 0, c2, 0, c3],
                  [c1, b1, c2, b2, c3, b3]]) / A / 2
    # print(B)
    if solve == 1:
        D = E / (1 - miu * miu) * np.array([[1, miu, 0],
                                            [miu, 1, 0],
                                            [0, 0, (1 - miu) / 2]])
    # elif solve == 2:
    #     D = E / (1 + miu) / (1 - 2 * miu) * np.array([[1 - miu, miu, 0],
    #                                                   [miu, 1 - miu, 0],
    #                                                   [0, 0, (1 - 2 * miu) / 2]])
    else:
        raise Exception("求解内容出错，请检查变量solve的值")
    S = np.dot(D, B)
    K_ele = t * A * np.dot(B.T, S)
    # print("单元刚度矩阵*16")
    # print(K_ele * 16)
    return K_ele, B, D


def assembly_triangle(K, K_ele, ele_node):
    # 点对点组装刚度矩阵
    # for i in range(len(ele_node)):
    #     for j in range(len(ele_node)):
    #         K[[ele_node[i] * 2, ele_node[i] * 2 + 1, ele_node[i] * 2, ele_node[i] * 2 + 1],
    #           [ele_node[j] * 2, ele_node[j] * 2, ele_node[j] * 2 + 1, ele_node[j] * 2 + 1]] \
    #             +=K_ele[[i * 2, i * 2 + 1, i * 2, i * 2 + 1], [j * 2, j * 2, j * 2 + 1, j * 2 + 1]]
    # print(K)
    # 向量化矩阵操作
    ele_node_ts = np.concatenate(([ele_node.reshape(len(ele_node), 1) * 2, ele_node.reshape(len(ele_node), 1) * 2 + 1]),
                                 axis=1).flatten()
    K[np.ix_(ele_node_ts, ele_node_ts)] += K_ele
    # K1[np.ix_(ele_node_ts, ele_node_ts)] += K_ele
    # print(np.array_equal(K, K1))
    return K


def outdata(u, Bs, Ds, ele_sum, ele_index, solve):
    stress = np.zeros([ele_sum, 3])
    strain = np.zeros([ele_sum, 6])
    node_stress = np.zeros([int(len(u.T) / 2), 3])
    # for ele in range(ele_sum):
    #     if solve == 1:
    #         stress[ele] = np.dot(np.dot(Ds[ele], Bs[ele]),
    #                              np.reshape(u.reshape(int(len(u.T) / 2), 2)[ele_index[ele], ...], [6, 1])).T
    #         node_stress[ele_index[ele], ...] += stress[ele, ...]
    #     # elif solve==2:
    #     #     strain[ele]= np.dot(Bs[ele],np.reshape(u.reshape(int(len(u.T)/2), 2)[ele_index[ele], ...],[6,1])).T
    #     else:
    #         raise Exception("求解内容出错，请检查变量solve的值")
    # 向量化部分循环操作
    if solve == 1:
        stress = np.matmul(np.matmul(Ds, Bs),
                           np.reshape(u.reshape(int(len(u.T) / 2), 2)[ele_index, ...],
                                      [ele_sum, 6, 1])).squeeze()
        for ele in range(ele_sum):
            node_stress[ele_index, ...] += stress[ele, ...]
    elif solve == 2:
        pass
    else:
        raise Exception("求解内容出错，请检查变量solve的值")
    unique_elements, counts = np.unique(ele_index, return_counts=True)
    node_stress[unique_elements, ...] /= (counts.reshape(len(counts), 1) * np.ones([1, 3]))
    print("单元应力")
    print(stress)
    print("节点应力")
    print(node_stress)
    # print(strain)
    return node_stress


def FEM(E, miu, t, node_ele, ele_index, f, BCs, solve):
    # 2D线性三角形单元
    node_sum = np.size(node_ele, 0)
    ele_sum = np.size(ele_index, 0)
    print("节点数：{}    单元数：{}".format(node_sum, ele_sum))
    print("节点坐标")
    print(node_ele.reshape(node_sum, 2))
    # if (ele_index.min() == 1) & (ele_index.max() == node_sum):
    #     ele_index = ele_index - 1
    #     ele_index=ele_index.astype(int)
    #     print("单元节点索引处理")
    # 检查节点索引是否从零开始，如果不是，索引自动减1
    if ele_index.min() == 1:
        ele_index = ele_index - 1
        ele_index = ele_index.astype(int)
        print("单元节点最小索引为1，自动处理")
    K = np.zeros([node_sum * 2, node_sum * 2])
    Bs = np.zeros([ele_sum, 3, 6])
    Ds = np.zeros([ele_sum, 3, 3])
    for ele in range(ele_sum):
        K_ele, Bs[ele], Ds[ele] = element_triangle(E, miu, t, node_ele, ele_index[ele], solve)
        K = assembly_triangle(K, K_ele, ele_index[ele])
    # print("刚度矩阵*16")
    # print(K * 16)
    F = np.zeros([1, node_sum * 2])
    for ele_f in f:
        F[[0, 0], [ele_f[0] * 2, ele_f[0] * 2 + 1]] = ele_f[1:]
    # print("载荷")
    # print(F)
    uf = np.zeros([1, node_sum * 2])
    u = np.zeros([1, node_sum * 2])
    for BC in BCs:
        uf[[0, 0], [BC[0] * 2, BC[0] * 2 + 1]] = BC[1:]
    # 行列删除法求约化后的刚度矩阵
    Kff = np.reshape(K[np.dot(uf.T == 0, uf == 0)], [np.sum(uf == 0), np.sum(uf == 0)])
    # 使用 solve 函数求解线性方程组
    u[uf == 0] = np.linalg.solve(Kff, F[uf == 0])
    print("位移")
    print(u.reshape(node_sum, 2))
    node_stress = outdata(u, Bs, Ds, ele_sum, ele_index, solve)
    print("有限元计算完成")
    return u.reshape(node_sum, 2), node_stress


if __name__ == '__main__':
    # 1.材料参数
    E = 1  # 弹性模量
    miu = 0  # 泊松比
    # 2.几何信息
    t = 1  # 厚度
    node_ele = np.array([[0, 0], [0, -2], [4, -1], [4, 0]])  # 节点坐标，按节点号排列(节点号必须为整数)
    ele_index = np.array([[1, 2, 0], [2, 3, 0]])  # 单元节点号，单元内逆时针排列
    # 3.荷载和边界条件
    f = np.array([[3, 0, -1]])  # 输入荷载信息，第一列为荷载作用的节点号，第二列为x方向荷载，第三列为y方向荷载
    BCs = np.array([[0, 1, 1], [1, 1, 1]])  # 约束信息，第一列为节点号，二~三列分别为x和y向的约束情况，1为固定，0为自由
    solve = 1  # 求解内容，1为应力
    FEM(E, miu, t, node_ele, ele_index, f, BCs, solve)
