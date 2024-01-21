"""
姓名: 王子聪
学号: S230200200

使用基于泊松方程的线弹性方程对跨度较短的桥梁进行有限元分析，
其中桥两侧受到固定支撑，桥面受到重力影响

桥的单元数为 5056个，在 GMSH 中画出

"""

import numpy as np

from wzc.utilities.plot import plot_results, plot_displacements
from wzc.mesh import Mesh
from wzc.fem.setup import FEMSetup


if __name__ == '__main__':
    # 读取网格
    mesh = Mesh('meshes/bridge.msh')
    # print(mesh.physical_groups_mapping)
    # 施加左右的约束
    mesh.set_boundary_condition(Mesh.BoundaryConditionType.DIRICHLET,
                                ['left_edge_bridge', 'right_edge_bridge'])
    # 添加重力约束
    mesh.set_boundary_condition(Mesh.BoundaryConditionType.NEUMANN,
                                ['up_left_bridge'])

    # mesh.draw()

    # 设置材料属性和重力大小
    fem = FEMSetup(
        mesh=mesh,
        rhs_func=lambda x: np.array([0, -9.81]),
        dirichlet_func=lambda x: np.array([0, 0]),
        neumann_func=lambda x: np.array([0, -9.81 * 1e2]),
        young_modulus=2.8e9,
        poisson_ratio=0.38,
    )
    # 求解
    results = fem.solve()

    # 计算形变和
    displacements = np.vstack((results[:mesh.nodes_num], results[mesh.nodes_num:])).T
    displacement_magnitudes = np.sqrt(displacements[:, 0] ** 2 + displacements[:, 1] ** 2)

    # 绘制结果图示
    plot_results(mesh, displacement_magnitudes)
    zoom_factor = 1e4
    plot_displacements(mesh, displacements * zoom_factor)
