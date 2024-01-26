"""
姓名：樊夫
学号：B2302S0112
程序目的：在梁的末端施加一个垂直向下的集中载荷2000 N，采用线性三角形单元计算梁的变形。
条件：假设材料是均质、各向同性的，且不考虑梁的自重，材料的 E=10e9 Pa nv=0.5
"""

import numpy as np

def element_stiffness_matrix(node_pos, E, miu):
    # 计算线性三角形单元的局部刚度矩阵
    B = np.zeros((3, 6))
    A = np.linalg.det([[1, node_pos[0, 0], node_pos[0, 1]],
                       [1, node_pos[1, 0], node_pos[1, 1]],
                       [1, node_pos[2, 0], node_pos[2, 1]]]) / 2
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        B[0, 2*i] = node_pos[j, 1] - node_pos[k, 1]
        B[1, 2*i+1] = node_pos[k, 0] - node_pos[j, 0]
        B[2, 2*i] = B[1, 2*i+1]
        B[2, 2*i+1] = B[0, 2*i]

    D = E / (1 - miu ** 2) * np.array([[1, miu, 0], [miu, 1, 0], [0, 0, (1 - miu) / 2]])

    return A * B.T @ D @ B

def assemble_stiffness_matrix(nodes, elements, E, miu):
    # 组装全局刚度矩阵
    N = np.zeros((2 * num_nodes, 2 * num_nodes))
    for element in elements:
        node_piece = element
        node_pos = nodes[node_piece]
        Ke = element_stiffness_matrix(node_pos, E, miu)
        for i in range(3):
            for j in range(3):
                N[2*node_piece[i]:2*node_piece[i]+2, 2*node_piece[j]:2*node_piece[j]+2] += Ke[2*i:2*i+2, 2*j:2*j+2]
    return N

def apply_boundary_conditions(N, nodes, fixed_node_indices):
    # 应用边界条件
    for node_idx in fixed_node_indices:
        for i in range(2):
            dof = 2 * node_idx + i
            N[dof, :] = 0
            N[:, dof] = 0
            N[dof, dof] = 1
    return N

def apply_loads(num_nodes, load_node_index, load_value):
    # 施加载荷
    F = np.zeros(2 * num_nodes)
    F[2 * load_node_index + 1] = -load_value  # 垂直向下的载荷
    return F

def solve_dis(K, F):
    # 求解节点位移
    dis = np.linalg.solve(K, F)
    return dis

# 梁的几何和网格参数
bar_len = 10  # 长度
bar_h = 1   # 高度
tri_num = 100  # 三角形元素的数量

# 材料属性
E = 10e9  # 杨氏模量 (Pa)
miu = 0.5   # 泊松比

# 生成节点和元素
nodes = []
elements = []
for i in range(tri_num + 1):
    for j in range(2):
        nodes.append([i * bar_len / bar_h, j * bar_h])

nodes = np.array(nodes)
num_nodes = len(nodes)

# 每个元素由3个节点定义
for i in range(tri_num):
    elements.append([i * 2, i * 2 + 1, i * 2 + 2])
    elements.append([i * 2 + 1, i * 2 + 3, i * 2 + 2])

elements = np.array(elements)

# 边界条件
N = assemble_stiffness_matrix(nodes, elements, E, miu)
N = apply_boundary_conditions(N, nodes, [0, 1])  # 假设第一个节点固定
F = apply_loads(num_nodes, num_nodes - 1, 2000)  # 在最后一个节点上施加2000 N的载荷
dis = solve_dis(N, F)

# 输出结果
print("节点位移矩阵:", dis)