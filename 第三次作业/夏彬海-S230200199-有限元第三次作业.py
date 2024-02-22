"""
姓名：夏彬海
学号：S230200199
程序目的：假设材料是均质、各向同性的，忽略梁的自重，假设在梁的末端施加一个垂直向下的集中载荷。求梁的形变
采用线性三角形单元
边界条件：固定支撑：程序中，第一个节点（节点索引为0）和第二个节点（节点索引为1）被设定为固定支撑。这意味着这些节点在两个方向（水平和垂直）上的位移都被限制为零。
         集中载荷：在梁的末端（最后一个节点，节点索引为num_nodes - 1）上施加了一个垂直向下的集中载荷。

"""


import numpy as np

# 梁的几何和网格参数
length = 10  # 梁的长度
height = 1   # 梁的高度
num_elements = 100  # 三角形元素的数量

# 材料属性
E = 210e9  # 杨氏模量 (Pa)
nu = 0.3   # 泊松比

# 生成节点和元素
nodes = []
elements = []
for i in range(num_elements + 1):
    for j in range(2):
        nodes.append([i * length / num_elements, j * height])

nodes = np.array(nodes)
num_nodes = len(nodes)

# 每个元素由3个节点定义
for i in range(num_elements):
    elements.append([i * 2, i * 2 + 1, i * 2 + 2])
    elements.append([i * 2 + 1, i * 2 + 3, i * 2 + 2])

elements = np.array(elements)

def element_stiffness_matrix(node_coords, E, nu):
    # 计算线性三角形单元的局部刚度矩阵
    # 简化的平面应力条件
    B = np.zeros((3, 6))
    A = np.linalg.det([[1, node_coords[0, 0], node_coords[0, 1]],
                       [1, node_coords[1, 0], node_coords[1, 1]],
                       [1, node_coords[2, 0], node_coords[2, 1]]]) / 2
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        B[0, 2*i] = node_coords[j, 1] - node_coords[k, 1]
        B[1, 2*i+1] = node_coords[k, 0] - node_coords[j, 0]
        B[2, 2*i] = B[1, 2*i+1]
        B[2, 2*i+1] = B[0, 2*i]

    D = E / (1 - nu ** 2) * np.array([[1, nu, 0],
                                      [nu, 1, 0],
                                      [0, 0, (1 - nu) / 2]])

    return A * B.T @ D @ B

def assemble_stiffness_matrix(nodes, elements, E, nu):
    # 组装全局刚度矩阵
    K = np.zeros((2 * num_nodes, 2 * num_nodes))
    for element in elements:
        node_indices = element
        node_coords = nodes[node_indices]
        Ke = element_stiffness_matrix(node_coords, E, nu)
        for i in range(3):
            for j in range(3):
                K[2*node_indices[i]:2*node_indices[i]+2, 2*node_indices[j]:2*node_indices[j]+2] += Ke[2*i:2*i+2, 2*j:2*j+2]
    return K

def apply_boundary_conditions(K, nodes, fixed_node_indices):
    # 应用边界条件
    for node_idx in fixed_node_indices:
        for i in range(2):
            dof = 2 * node_idx + i
            K[dof, :] = 0
            K[:, dof] = 0
            K[dof, dof] = 1
    return K

def apply_loads(num_nodes, load_node_index, load_value):
    # 施加载荷
    F = np.zeros(2 * num_nodes)
    F[2 * load_node_index + 1] = -load_value  # 垂直向下的载荷
    return F

def solve_displacements(K, F):
    # 求解节点位移
    displacements = np.linalg.solve(K, F)
    return displacements

# 主要的FEA过程
K = assemble_stiffness_matrix(nodes, elements, E, nu)
K = apply_boundary_conditions(K, nodes, [0, 1])  # 假设第一个节点固定
F = apply_loads(num_nodes, num_nodes - 1, 1000)  # 在最后一个节点上施加1000 N的载荷
displacements = solve_displacements(K, F)

# 输出结果
print("节点位移:", displacements)
