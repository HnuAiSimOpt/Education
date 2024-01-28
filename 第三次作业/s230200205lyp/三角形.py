"""
问题描述:单一力下的节点位移、单元应变和单元应力问题
方法：采用线性三角单元
姓名：李炎朴S230200205
"""

import numpy as np
import scipy as s


def create_K(EL, PO, TK, NODE_SUM, ELE_SUM, NODE, ELE):  # 生成整体刚度矩阵
    ASK = np.zeros((2 * NODE_SUM, 2 * NODE_SUM))
    for n in range(1, ELE_SUM + 1):
        x = np.array([0] * 6)
        y = np.array([0] * 6)
        for i in range(3):
            x[i] = NODE[ELE[n - 1, i] - 1, 0]
            y[i] = NODE[ELE[n - 1, i] - 1, 1]
        x[3:6] = x[0:3]
        y[3:6] = y[0:3]

        MA = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])
        A = 0.5 * np.linalg.det(MA)

        a = np.array([0] * 3)
        b = np.array([0] * 3)
        c = np.array([0] * 3)
        for i in range(3):
            a[i] = x[i + 1] * y[i + 2] - x[i + 2] * y[i + 1]
            b[i] = y[i + 1] - y[i + 2]
            c[i] = -x[i + 1] + x[i + 2]

        EK = np.zeros((6, 6))
        for i in range(3):
            for j in range(3):
                k1 = b[i] * b[j] + (1 - PO) / 2 * c[i] * c[j]
                k2 = PO * c[i] * b[j] + (1 - PO) / 2 * b[i] * c[j]
                k3 = PO * b[i] * c[j] + (1 - PO) / 2 * c[i] * b[j]
                k4 = c[i] * c[j] + (1 - PO) / 2 * b[i] * b[j]
                EK[2 * i:2 * i + 2, 2 * j:2 * j + 2] = EL * TK / (4 * (1 - PO ** 2) * A) * np.array(
                    [[k1, k3], [k2, k4]])

        G = np.zeros((6, 2 * NODE_SUM))
        for i in range(3):
            G[2 * i, 2 * ELE[n - 1, i] - 2] = 1
            G[2 * i + 1, 2 * ELE[n - 1, i] - 1] = 1
        GT = np.transpose(G)
        ask = np.dot(np.dot(GT, EK), G)
        ASK = ASK + ask
    return ASK


def create_P(P_inp, NODE_SUM):  # 由输入的节点荷载信息生成荷载向量P
    P = np.zeros((2 * NODE_SUM, 1))
    k = 1
    P2 = P_inp.astype(int)
    for i in P2[:, 0]:
        P[2 * i - 2] = P_inp[k - 1, 1]
        P[2 * i - 1] = P_inp[k - 1, 2]
        k += 1
    return P


def do_BC(BC):  # 处理边界条件，置1法对角元素
    global ASK
    global P
    k = 1
    for i in BC[:, 0]:
        if BC[k - 1, 1] == 1:
            P[2 * i - 2] = 0
        if BC[k - 1, 2] == 1:
            P[2 * i - 1] = 0
        k += 1
    k = 1
    for i in BC[:, 0]:
        if BC[k - 1, 1] == 1:
            ASK[2 * i - 2, :] = 0
            ASK[:, 2 * i - 2] = 0
            ASK[2 * i - 2, 2 * i - 2] = 1
        if BC[k - 1, 2] == 1:
            ASK[2 * i - 1, :] = 0
            ASK[:, 2 * i - 1] = 0
            ASK[2 * i - 1, 2 * i - 1] = 1
        k += 1
    return


def solve_Uout(U, NODE_SUM):
    Uout = np.zeros((NODE_SUM, 3))
    for i in range(NODE_SUM):
        Uout[i, 0] = i + 1
        Uout[i, 1] = U[2 * i]
        Uout[i, 2] = U[2 * i + 1]
    return Uout


def solve_B(U, NODE, ELE, ELE_SUM):
    B = np.zeros((ELE_SUM, 4))
    for n in range(1, ELE_SUM + 1):

        ue = np.zeros((6, 1))
        for i in range(3):
            ue[2 * i, 0] = U[2 * ELE[n - 1, i] - 2]
            ue[2 * i + 1, 0] = U[2 * ELE[n - 1, i] - 1]

        x = np.array([0] * 6)
        y = np.array([0] * 6)
        for i in range(3):
            x[i] = NODE[ELE[n - 1, i] - 1, 0]
            y[i] = NODE[ELE[n - 1, i] - 1, 1]
        x[3:6] = x[0:3]
        y[3:6] = y[0:3]

        a = np.array([0] * 3)
        b = np.array([0] * 3)
        c = np.array([0] * 3)
        for i in range(3):
            a[i] = x[i + 1] * y[i + 2] - x[i + 2] * y[i + 1]
            b[i] = y[i + 1] - y[i + 2]
            c[i] = -x[i + 1] + x[i + 2]

        MA = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])
        A = 0.5 * np.linalg.det(MA)

        BE = 1 / 2 * A * np.array(
            [[b[0], 0, b[1], 0, b[2], 0], [0, c[0], 0, c[1], 0, c[2]], [c[0], b[0], c[1], b[1], c[2], b[2]]])
        be = np.dot(BE, ue)

        B[n - 1, :] = [n, be[0, 0], be[1, 0], be[2, 0]]
    return B


def solve_S(U, NODE, ELE, ELE_SUM, EL, PO):
    S = np.zeros((ELE_SUM, 4))
    for n in range(1, ELE_SUM + 1):

        ue = np.zeros((6, 1))
        for i in range(3):
            ue[2 * i, 0] = U[2 * ELE[n - 1, i] - 2]
            ue[2 * i + 1, 0] = U[2 * ELE[n - 1, i] - 1]

        x = np.array([0] * 6)
        y = np.array([0] * 6)
        for i in range(3):
            x[i] = NODE[ELE[n - 1, i] - 1, 0]
            y[i] = NODE[ELE[n - 1, i] - 1, 1]
        x[3:6] = x[0:3]
        y[3:6] = y[0:3]

        a = np.array([0] * 3)
        b = np.array([0] * 3)
        c = np.array([0] * 3)
        for i in range(3):
            a[i] = x[i + 1] * y[i + 2] - x[i + 2] * y[i + 1]
            b[i] = y[i + 1] - y[i + 2]
            c[i] = -x[i + 1] + x[i + 2]

        MA = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])
        A = 0.5 * np.linalg.det(MA)

        SE = np.zeros((3, 6))
        for i in range(3):
            SE[:, 2 * i:2 * i + 2] = EL / (2 * (1 - PO ** 2) * A) * np.array(
                [[b[i], PO * c[i]], [PO * b[i], c[i]], [(1 - PO) / 2 * c[i], (1 - PO) / 2 * b[i]]])
        se = np.dot(SE, ue)

        S[n - 1, :] = [n, se[0, 0], se[1, 0], se[2, 0]]
    return S


# 1.材料参数
EL = 20  # 弹性模量
PO = 0.1  # 泊松比

# 2.几何信息
TK = 1  # 厚度
NODE_SUM = 5  # 节点数目
ELE_SUM = 3  # 单元数目
NODE = np.array([[0, 1], [0, 0], [2, 0], [2, 1], [4, 1]])  # 节点坐标
ELE = np.array([[1, 2, 3], [3, 4, 1], [3, 5, 4]])  # 单元节点号

# 3.荷载和边界条件
P_inp = np.array([[5, 0, -1]])  # 输入荷载
BC = np.array([[1, 1, 1], [2, 1, 1]])  # 约束

ASK = create_K(EL, PO, TK, NODE_SUM, ELE_SUM, NODE, ELE)  # 生成结构刚度矩阵
P = create_P(P_inp, NODE_SUM)  # 生成荷载矩阵
do_BC(BC)  # 处理边界条件
U = s.sparse.linalg.spsolve(ASK, P)  # 解大型稀疏矩阵方程

Uout = solve_Uout(U, NODE_SUM)
B = solve_B(U, NODE, ELE, ELE_SUM)  # 得到单元应变信息矩阵
S = solve_S(U, NODE, ELE, ELE_SUM, EL, PO)  # 得到单元应力信息矩阵

print('##节点位移##\n节点号\tx\ty')
print(Uout)
print('\n##单元应变##\n单元号\tεx\tεy\tγxy')
print(B)
print('\n##单元应力##\n单元号\tσx\tσy\tτxy')
print(S)
