
import numpy as np
from scipy.linalg import solve

import warnings
warnings.filterwarnings("ignore")

# 材料属性和几何参数
E = 210e9  # 弹性模量，单位Pa
t = 0.01   # 单元厚度，单位m
v = 0.3    # 泊松比
h = 0.1    # 梁高度，单位m
L = 1.0    # 梁长度，单位m

# 创建一个简单的梁，用两个三角形单元
nodes = np.array([
    [0, 0],  # 节点1
    [L, 0],  # 节点2
    [L, h],  # 节点3
    [0, h]   # 节点4
])

# 单元连接性
elements = np.array([
    [1, 2, 4],  # 单元1
    [2, 3, 4]   # 单元2
]) - 1

# 边界条件
# 固定节点1和4的所有自由度
fixed_dofs = [0, 1, 6, 7]

# 载荷
# 假设在节点2和3上施加一个向下的集中载荷
P = -10000  # 集中载荷大小，单位N
force = np.zeros((len(nodes) * 2, 1))
force[2] = P/2  # 节点2的y方向
force[4] = P/2  # 节点3的y方向

# 初始化全局刚度矩阵和全局载荷
K_global = np.zeros((len(nodes) * 2, len(nodes) * 2))
F_global = np.zeros((len(nodes) * 2, 1))

# 计算单元刚度矩阵和组装全局刚度矩阵
for el in elements:
    # 提取单元节点坐标
    coords = nodes[el, :]
    # 计算单元面积
    A = 0.5 * np.linalg.det(np.array([[1, coords[0, 0], coords[0, 1]],
                                       [1, coords[1, 0], coords[1, 1]],
                                       [1, coords[2, 0], coords[2, 1]]]))
    # 计算单元刚度矩阵
    B = np.array([[coords[1, 1] - coords[2, 1], 0, coords[2, 1] - coords[0, 1], 0, coords[0, 1] - coords[1, 1], 0],
                  [0, coords[2, 0] - coords[1, 0], 0, coords[0, 0] - coords[2, 0], 0, coords[1, 0] - coords[0, 0]],
                  [coords[2, 0] - coords[1, 0], coords[1, 1] - coords[2, 1], coords[0, 0] - coords[2, 0], coords[2, 1] - coords[0, 1], coords[1, 0] - coords[0, 0], coords[0, 1] - coords[1, 1]]]) / (2 * A)
    D = E / (1 - v**2) * np.array([[1, v, 0],
                                   [v, 1, 0],
                                   [0, 0, (1 - v) / 2]])
    k_el = t * A * np.dot(np.dot(B.T, D), B)

    # 组装全局刚度矩阵
    dofs = np.array([2 * el[0], 2 * el[0] + 1, 2 * el[1], 2 * el[1] + 1, 2 * el[2], 2 * el[2] + 1])
    for i in range(6):
        for j in range(6):
            K_global[dofs[i], dofs[j]] += k_el[i, j]

# 应用边界条件
K_reduced = np.delete(np.delete(K_global, fixed_dofs, axis=0), fixed_dofs, axis=1)
F_reduced = np.delete(force, fixed_dofs, axis=0)

# 解线性方程组得到位移
U_reduced = solve(K_reduced, F_reduced)

# 插入固定自由度的位移（
U = np.zeros((len(nodes) * 2, 1))
U[np.setdiff1d(range(len(nodes) * 2), fixed_dofs)] = U_reduced

# 计算应力和应变
stress = np.zeros((len(elements), 3))  # 每个单元3个应力分量
strain = np.zeros((len(elements), 3))  # 每个单元3个应变分量
for i, el in enumerate(elements):
    dofs = np.array([2 * el[0], 2 * el[0] + 1, 2 * el[1], 2 * el[1] + 1, 2 * el[2], 2 * el[2] + 1])
    U_el = U[dofs, :]
    coords = nodes[el, :]
    A = 0.5 * np.linalg.det(np.array([[1, coords[0, 0], coords[0, 1]],
                                       [1, coords[1, 0], coords[1, 1]],
                                       [1, coords[2, 0], coords[2, 1]]]))
    B = np.array([[coords[1, 1] - coords[2, 1], 0, coords[2, 1] - coords[0, 1], 0, coords[0, 1] - coords[1, 1], 0],
                  [0, coords[2, 0] - coords[1, 0], 0, coords[0, 0] - coords[2, 0], 0, coords[1, 0] - coords[0, 0]],
                  [coords[2, 0] - coords[1, 0], coords[1, 1] - coords[2, 1], coords[0, 0] - coords[2, 0], coords[2, 1] - coords[0, 1], coords[1, 0] - coords[0, 0], coords[0, 1] - coords[1, 1]]]) / (2 * A)
    strain[i, :] = np.dot(B, U_el).flatten()
    stress[i, :] = np.dot(D, strain[i, :])

# 输出结果
print("节点位移:")
print(U)
print("单元应力:")
print(stress)
print("单元应变:")
print(strain)








