import pygmsh
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.colors as mcolors

plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号

import Isoparametric_quadrangular_elelments_S230200206 as iqe
import linear_triangular_elements_S200600206 as lte


def quads_to_tris(quads):
    # 四边形拆为三角形
    # tris = [[None for j in range(3)] for i in range(2 * len(quads))]
    # for i in range(len(quads)):
    #     j = 2 * i
    #     n0 = quads[i][0]
    #     n1 = quads[i][1]
    #     n2 = quads[i][2]
    #     n3 = quads[i][3]
    #     tris[j][0] = n0
    #     tris[j][1] = n1
    #     tris[j][2] = n2
    #     tris[j + 1][0] = n2
    #     tris[j + 1][1] = n3
    #     tris[j + 1][2] = n0
    # 向量化循环操作
    tris = np.empty([quads.shape[0] * 2, 3], dtype=quads.dtype)
    tris[::2, :] = quads[:, [0, 1, 2]]
    tris[1::2, :] = quads[:, [2, 3, 0]]
    return np.array(tris)


def plot_fem_mesh(nodes_x, nodes_y, elements, c):
    for element in elements:
        x = [nodes_x[element[i]] for i in range(len(element))]
        y = [nodes_y[element[i]] for i in range(len(element))]
        plt.fill(x, y, edgecolor=c, fill=False)


with pygmsh.occ.Geometry() as geom:
    geom.characteristic_length_min = 5
    geom.characteristic_length_max = 5
    # 定义一个矩形
    rectangle = geom.add_rectangle([0.0, 0.0, 0.0], 200.0, 100.0)

    # 定义梯形
    # geom.add_polygon(
    #     [
    #         [0, 0],
    #         [0, -2],
    #         [4, -1],
    #         [4, 0]
    #     ],
    #     mesh_size=5,
    # )
    # 切割出小孔
    disk1 = geom.add_disk([100.0, 50.0, 0.0], 10)
    flat = geom.boolean_difference(rectangle, disk1)
    # 三维拉伸
    # geom.extrude(flat, [0.0, 0.0, 1])
    # 生成网格,11为Frontal-Delaunay for Quad
    # mesh = geom.generate_mesh(2, algorithm=11)
    mesh = geom.generate_mesh(2, algorithm=8)
    # 1.材料参数
    E = 200 * 1000  # 弹性模量
    miu = 0.3  # 泊松比
    # 2.几何信息
    t = 1  # 厚度
    node_ele = np.array(mesh.points)[:, 0:2]
    ele_index = np.array(mesh.cells[1].data)
    # 3.荷载和边界条件
    f_i = np.where((node_ele[:, 0] == 200))[0][:, np.newaxis]
    f = np.reshape(np.concatenate((f_i, np.tile(np.array([2000, 0]), (f_i.shape[0], 1))), axis=1), [f_i.shape[0], 3])
    BCs_i = np.array(np.where(node_ele[:, 0] == 0)).T
    BCs = np.ones([BCs_i.size, 3])
    BCs[:, 0] = BCs_i.reshape([BCs_i.size, ])
    BCs = BCs.astype(int)
    solve = 1  # 求解内容，1为应力
    gauss = 2  # 高斯积分形式，2为全积分

    # 绘制节点的位移
    plt.figure()
    # 创建并添加节点
    if ele_index.shape[1] == 4:
        triang = mtri.Triangulation(node_ele[:, 0], node_ele[:, 1], quads_to_tris(ele_index))
    elif ele_index.shape[1] == 3:
        triang = mtri.Triangulation(node_ele[:, 0], node_ele[:, 1], ele_index)
    # plt.triplot(triang, 'go-', lw=1)
    # plt.triplot(triang2, 'ro-', lw=1)
    plot_fem_mesh(node_ele[:, 0], node_ele[:, 1], ele_index, 'black')
    plt.show()
    if ele_index.shape[1] == 4:
        u, node_stress = iqe.FEM(E, miu, t, node_ele, ele_index, f, BCs, gauss, solve)
    elif ele_index.shape[1] == 3:
        u, node_stress = lte.FEM(E, miu, t, node_ele, ele_index, f, BCs, solve)

    node_sum = np.size(node_ele, 0)
    node_ele_u = node_ele + u
    stress = np.linalg.norm(node_stress, axis=1)
    # 可视化
    print("结果可视化")

    # 绘制节点的位移
    plt.figure()
    # 创建并添加节点
    if ele_index.shape[1] == 4:
        triang = mtri.Triangulation(node_ele[:, 0], node_ele[:, 1], quads_to_tris(ele_index))
        triang2 = mtri.Triangulation(node_ele_u[:, 0], node_ele_u[:, 1], quads_to_tris(ele_index))
    elif ele_index.shape[1] == 3:
        triang = mtri.Triangulation(node_ele[:, 0], node_ele[:, 1], ele_index)
        triang2 = mtri.Triangulation(node_ele_u[:, 0], node_ele_u[:, 1], ele_index)
    # plt.triplot(triang, 'go-', lw=1)
    # plt.triplot(triang2, 'ro-', lw=1)
    plot_fem_mesh(node_ele[:, 0], node_ele[:, 1], ele_index, 'black')
    plot_fem_mesh(node_ele_u[:, 0], node_ele_u[:, 1], ele_index, 'red')
    if ele_index.shape[1] == 4:
        plt.title('四边形-位移图')
    elif ele_index.shape[1] == 3:
        plt.title('三角形-位移图')
    plt.show()
    # 节点应力图
    plt.figure()
    fig = plt.tricontourf(triang2, stress, cmap="jet")
    plt.colorbar(fig, location='right', orientation='vertical')
    if ele_index.shape[1] == 4:
        plt.title('四边形-应力图')
    elif ele_index.shape[1] == 3:
        plt.title('三角形-应力图')
    plt.show()
