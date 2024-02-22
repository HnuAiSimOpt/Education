import numpy as np

'''
学号：S230200160
姓名：陈桓隆
题目：板长为2m，宽为0.1m，上端面受均布载荷q=1MPa，两端面固定，弹性模量为2.1e4Pa，泊松比为0.3，板厚为0.01m，计算其应力和变形
'''
# 板的几何和网格参数
length = 2.0  # 板的长度
width = 0.1  # 板的宽度
thickness = 0.01  # 板的厚度
num_elements = 200  # 三角形元素的数量

# 材料属性
E = 2.1e10  # 弹性模量 (Pa)
nu = 0.3  # 泊松比

# 载荷属性
q = 1e6  # 均布载荷 (Pa)

# 生成节点和元素
nodes = []
elements = []
for i in range(num_elements + 1):
    for j in range(2):
        nodes.append([i * length / num_elements, j * width])

nodes = np.array(nodes)
num_nodes = len(nodes)

# 每个元素由3个节点定义
for i in range(num_elements):
    elements.append([i * 2, i * 2 + 1, i * 2 + 2])
    elements.append([i * 2 + 1, i * 2 + 3, i * 2 + 2])

elements = np.array(elements)


def element_stiffness_matrix(node_coords, E, nu, thickness):
    # 计算线性三角形单元的局部刚度矩阵
    B = np.zeros((3, 6))
    A = np.linalg.det([[1, node_coords[0, 0], node_coords[0, 1]],
                       [1, node_coords[1, 0], node_coords[1, 1]],
                       [1, node_coords[2, 0], node_coords[2, 1]]]) / 2
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        B[0, 2 * i] = node_coords[j, 1] - node_coords[k, 1]
        B[1, 2 * i + 1] = node_coords[k, 0] - node_coords[j, 0]
        B[2, 2 * i] = B[1, 2 * i + 1]
        B[2, 2 * i + 1] = B[0, 2 * i]

    D = E / (1 - nu ** 2) * np.array([[1, nu, 0],
                                      [nu, 1, 0],
                                      [0, 0, (1 - nu) / 2]])

    return thickness * A * B.T @ D @ B


def assemble_stiffness_matrix(nodes, elements, E, nu, thickness):
    # 组装全局刚度矩阵
    K = np.zeros((2 * num_nodes, 2 * num_nodes))
    for element in elements:
        node_indices = element
        node_coords = nodes[node_indices]
        Ke = element_stiffness_matrix(node_coords, E, nu, thickness)
        for i in range(3):
            for j in range(3):
                K[2 * node_indices[i]:2 * node_indices[i] + 2, 2 * node_indices[j]:2 * node_indices[j] + 2] += Ke[
                                                                                                               2 * i:2 * i + 2,
                                                                                                               2 * j:2 * j + 2]
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


def apply_loads(num_nodes, q, thickness):
    # 施加载荷
    F = np.zeros(2 * num_nodes)
    for i in range(num_nodes):
        F[2 * i + 1] = -q * length * thickness / num_elements  # 在每个节点上施加均布载荷
    return F


def solve_displacements(K, F):
    # 求解节点位移
    displacements = np.linalg.solve(K, F)
    return displacements


# 主要的FEA过程
K = assemble_stiffness_matrix(nodes, elements, E, nu, thickness)
K = apply_boundary_conditions(K, nodes, [0, num_elements])  # 两端面固定
F = apply_loads(num_nodes, q, thickness)  # 上端面施加均布载荷
displacements = solve_displacements(K, F)

# 输出结果
print("节点位移:", displacements)
