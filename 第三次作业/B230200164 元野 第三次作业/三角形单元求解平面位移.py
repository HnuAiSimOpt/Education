#
import numpy as np
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import spsolve

# B230200164 元野
# 程序目的：利用三角形单元解决平面位移问题
# 单元类型：线性三角形单元
# 边界条件：边界条件固定左边和底部界面，上界面受到均匀的向下压力，并且假设右边界自由。

# 设置材料属性，这里我们使用钢铁的典型属性
E = 210e9  # 杨氏模量，单位：Pa
nu = 0.3   # 泊松比
thickness = 0.01  # 假设材料厚度1cm，对于平面应变来说，这个值仅影响力的大小，不影响分布

# 考虑的矩形区域
width = 2.0
height = 1.0

# 这里我们将矩形区域划分为若干个三角形元素，这部分通常由网格生成器完成
# 暂时我们假设网格由4个节点和2个元素组成
nodes = np.array([
    [0, 0], 
    [width, 0], 
    [width, height], 
    [0, height]
])

elements = np.array([
    [0, 1, 3], 
    [1, 2, 3]
])

# 每个节点的自由度（这里是二维问题）
dof_per_node = 2

# 计算平面应变状态的刚度矩阵
D = E / (1 - nu**2) * np.array([
    [1, nu, 0],
    [nu, 1, 0],
    [0, 0, (1 - nu) / 2]
])

# 形函数的局部导数
def dN_local():
    return np.array([
        [-1, 1, 0],
        [-1, 0, 1]
    ])

# 组装全局刚度矩阵
def assemble_global_stiffness(N, E, num_dof):
    K_global = lil_matrix((num_dof, num_dof))  # 使用LIL格式，因为我们将对刚度矩阵进行频繁的更改
    for el in elements:
        coord = nodes[el, :]
        K_element = element_stiffness_matrix(D, coord, thickness)
        assembly_element(K_global, K_element, el)
    return csc_matrix(K_global)  # 转换为CSC格式以进行求解

# 计算雅可比矩阵
def jacobian(coord):
    dN = dN_local()
    J = np.dot(dN, coord)
    return J

# 计算雅可比矩阵的逆和行列式
def inv_jacobian(coord):
    J = jacobian(coord)
    detJ = np.linalg.det(J)
    invJ = np.linalg.inv(J)
    return invJ, detJ

# 计算B矩阵
def B_matrix(coord):
    invJ, detJ = inv_jacobian(coord)
    dN = dN_local()
    dN_dx_dy = np.dot(invJ, dN)  # 在全局坐标中的导数
    B = np.zeros((3, 6))
    # 填充B矩阵
    for i in range(3):
        B[0, i*2] = dN_dx_dy[0, i]
        B[1, i*2+1] = dN_dx_dy[1, i]
        B[2, i*2] = dN_dx_dy[1, i]
        B[2, i*2+1] = dN_dx_dy[0, i]
    return B, detJ
    
# 完成原来省略掉的单元刚度矩阵计算函数
def element_stiffness_matrix(D, coord, thickness):
    B, detJ = B_matrix(coord)
    K_element = thickness * detJ * np.dot(np.dot(B.T, D), B)
    return K_element

# 元素刚度矩阵装配到全局刚度矩阵
def assembly_element(K_global, K_element, element_nodes):
    for i in range(len(element_nodes)):
        for j in range(len(element_nodes)):
            node_i = element_nodes[i]
            node_j = element_nodes[j]
            K_global[node_i*2:node_i*2+2, node_j*2:node_j*2+2] += K_element[i*2:i*2+2, j*2:j*2+2]

# 应用力和边界条件
def apply_bc_force(K_global, F_global):
    # 应用力
    # 假设上边界受到沿着y轴方向的均布荷载，力的总量为10000N（再次假设）
    F_global[-2] = -10000 / 2  # 分配到右上角节点
    F_global[-4] = -10000 / 2  # 分配到左上角节点

    # 应用边界条件
    for node_id in [0, 1]:  # 假设左下角节点固定不动
        for dof in [0, 1]:  # 两个自由度，X和Y
            dof_id = node_id * dof_per_node + dof
            K_global[dof_id, :] = 0
            K_global[:, dof_id] = 0
            K_global[dof_id, dof_id] = 1  # 给对角线元素赋值为1，这样可以确保矩阵求解时不会出现问题
            F_global[dof_id] = 0  # 边界上的力设为0

    return K_global, F_global

# 主程序
num_nodes = nodes.shape[0]
num_dof = num_nodes * dof_per_node
K_global = assemble_global_stiffness(nodes, elements, num_dof)
F_global = np.zeros(num_dof)

# 应用力和边界条件
K_global, F_global = apply_bc_force(K_global, F_global)

# 求解位移
displacements = spsolve(K_global.tocsr(), F_global)

print("Displacements:")
print(displacements)


