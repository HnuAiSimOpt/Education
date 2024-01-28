""""
线性三角形有限元计算算法

姓名 王子聪
学号 S230200200

本算法用于解决任意正四边形的有限元计算，
计算时需要指定四边形的长和宽，
材料的 E 和 nv 参数，
以及约束和附加的外力

本程序使用线性三角形法计算，
因此需要提前划分三角形网格

本算法的测试用例是一个 1X2 的正长方形，
其中宽边的一侧完全固定，
在另一侧的中点处施加一个 1N 的垂直外力，
材料的 E=10 nv=0.5
"""

import numpy as np
# 画图使用，不参与计算
import matplotlib.pyplot as plt
import matplotlib.tri as tri


class Material:
    def __init__(self, e: int | float, nv: int | float):
        """
        根据材料的 E 和 nv 计算本构矩阵 D
        """
        self.E = e
        self.nv = nv
        self.D = self.get_d()
        return

    def get_d(self) -> np.ndarray:
        """
        :return:
            返回材料的本构矩阵
        """
        return self.E / (1.0 - self.nv ** 2) * np.array([
            [1.0, self.nv, 0.0],
            [self.nv, 1.0, 0.0],
            [0.0, 0.0, 0.5 * (1.0 - self.nv)],
        ])


class TriangleElement:
    def __init__(self, material: Material, positions: np.ndarray):
        """
        根据指定的材料和三角的顶点坐标，
        计算 B 矩阵和刚度矩阵

        :param material:
            材料属性
        :param positions:
            三角形的顶点坐标
        """
        self.x = positions.T[0]
        self.y = positions.T[1]
        self.node_set = np.array([0, 1, 2], dtype=np.uint64)  #: 初始的节点信息
        self.D = material.D

        self.area_2 = self.get_area_2()
        self.B = self.get_b()
        self.K_element = self.element_integrate()
        return

    def update_node_set(self, node_set: np.ndarray):
        """
        更新节点数据
        """
        self.node_set = node_set
        return

    def get_area_2(self) -> np.ndarray:
        """
        :return:
            两倍三角形面积
        """
        matrix = np.array([
            [1.0, self.x[0], self.y[0]],
            [1.0, self.x[1], self.y[1]],
            [1.0, self.x[2], self.y[2]],
        ])
        return abs(np.linalg.det(matrix))

    def get_b(self) -> np.ndarray:
        """
        :return:
            B 矩阵
        """
        return 0.5 / self.area_2 * np.array([
            [self.y[1] - self.y[2], 0.0, self.y[2] - self.y[0],
             0.0, self.y[0] - self.y[1], 0.0],
            [0.0, self.x[2] - self.x[1], 0.0,
             self.x[0] - self.x[2], 0.0, self.x[1] - self.x[0]],
            [self.x[2] - self.x[1], self.y[1] - self.y[2], self.x[0] - self.x[2],
             self.y[2] - self.y[0], self.x[1] - self.x[0], self.y[0] - self.y[1]],
        ])

    def element_integrate(self) -> np.ndarray:
        """
        :return:
            拼接刚度矩阵
        """
        return self.area_2 * self.B.T @ self.D @ self.B

    def get_strain(self, deform_local) -> np.ndarray:
        """
        :return:
            三角形单元的压力值
        """
        return self.B @ deform_local


class Fea:
    def __init__(self, fea_mesh):
        """
        :param fea_mesh:
            三角形网格
        """
        self.mesh = fea_mesh

        # 处理三角形网格信息
        self.nodes = fea_mesh.points[:, :2]
        self.elements = fea_mesh.cells_dict['triangle']
        self.materials = fea_mesh.cell_data

        # 初始化变量
        self.len_reduce = None
        self.deform_free_index = None
        self.deform_reduce = None
        self.K_reduce = None
        self.force_reduce = None
        self.stress_list = None
        self.strain_list = None
        self.len_global = 2 * len(self.nodes)
        self.deform = np.zeros(self.len_global)
        self.force = np.zeros(self.len_global)

        # 初始化总体刚度矩阵
        self.K = None
        self.get_k()

        # 初始化边界信息
        self.x_fix = {}
        self.y_fix = {}
        # 初始化载荷信息
        self.f_given = {}

        # 初始化有限元计算结果字典
        self.result_dict = None
        self.init_show_dict()

        return

    def init_show_dict(self):
        """
        初始化结果字典
        """
        self.result_dict = {'position': {'x': [], 'y': []},
                            'deform': {'Ux': [], 'Uy': []},
                            'force': {'Fx': [], 'Fy': []},
                            'strain': {'e11': [], 'e22': [], 'e12': []},
                            'stress': {'S11': [], 'S22': [], 'S12': []},
                            }
        self.result_dict['position']['x'] = self.nodes.T[0]
        self.result_dict['position']['y'] = self.nodes.T[1]
        return

    @staticmethod
    def set_material(fea_mesh, material: Material):
        """
        设置三角形单元的材料属性

        :param fea_mesh:
            网格信息
        :param material:
            材料属性
        """
        for i, _ in enumerate(mesh.cells_dict['triangle']):
            fea_mesh.cell_data.update({i: material})
        return

    def set_conditions(self, x_fix: dict, y_fix: dict, f_given: dict):
        """
        设置边界约束和外部载荷

        :param x_fix:
            X 方向上的固定载荷的节点编号
        :param y_fix:
            Y 方向上的固定载荷的节点编号
        :param f_given:
            外部载荷作用的节点编号以及力沿坐标轴的分量
        """
        self.x_fix.update(x_fix)
        self.y_fix.update(y_fix)
        self.f_given.update(f_given)
        return

    def instant_element(self, element_index: int, node_set: np.ndarray) -> TriangleElement:
        """
        根据节点信息和索引坐标填充线性三角形单元

        :param element_index:
            索引坐标
        :param node_set:
            节点信息

        :return:
            填充完的线性三角形单元
        """
        material = self.materials[element_index]
        positions = np.array([self.nodes[node_set[j]][:] for j in range(3)])
        return TriangleElement(material, positions)

    def get_ke2k(self, element: TriangleElement):
        """
        填充单元刚度矩阵

        :param element:
            三角性单元
        """
        k_element = element.K_element
        node_set = element.node_set
        deform_global_index = np.array([[2 * node_set[i], 2 * node_set[i] + 1] for i in range(3)],
                                       dtype=np.uint64).reshape(-1)
        for i_local, i_global in enumerate(deform_global_index):
            for j_local, j_global in enumerate(deform_global_index):
                self.K[i_global, j_global] += k_element[i_local, j_local]
        return

    def get_k(self):
        """
        填充整体刚度矩阵
        """
        self.K = np.zeros((self.len_global, self.len_global))
        for element_index, node_set in enumerate(self.elements):
            element = self.instant_element(element_index, node_set)
            element.update_node_set(node_set)
            self.get_ke2k(element)
        return

    def get_fem(self):
        """
        根据约束和载荷信息，
        以及材料属性矩阵，
        生成待求解的 FEM 方程
        """
        self.deform_free_index = []
        for node, _ in enumerate(self.nodes):
            if node not in self.x_fix:
                self.deform_free_index.append(2 * node)
            if node not in self.y_fix:
                self.deform_free_index.append(2 * node + 1)
        self.len_reduce = len(self.deform_free_index)

        for node in self.x_fix:
            self.deform[2 * node] = self.x_fix[node]
        for node in self.y_fix:
            self.deform[2 * node + 1] = self.y_fix[node]
        for node in self.f_given:
            self.force[2 * node] = self.f_given[node][0]
            self.force[2 * node + 1] = self.f_given[node][1]

        self.deform_reduce = np.zeros(self.len_reduce)
        self.force_reduce = np.zeros(self.len_reduce)
        self.K_reduce = np.zeros((self.len_reduce, self.len_reduce))
        for i_reduce, i_global in enumerate(self.deform_free_index):
            self.force_reduce[i_reduce] = self.force[i_global]
        for i_reduce, i_global in enumerate(self.deform_free_index):
            for j_reduce, j_global in enumerate(self.deform_free_index):
                self.K_reduce[i_reduce, j_reduce] = self.K[i_global, j_global]
        return

    def solve_fem(self):
        """
        使用 numpy 对 FEM 方程进行求解
        """
        self.deform_reduce = np.linalg.solve(self.K_reduce, self.force_reduce)
        return

    def update_global_variables(self):
        """
        更新变量
        """
        for i_reduce, i_global in enumerate(self.deform_free_index):
            self.deform[i_global] = self.deform_reduce[i_reduce]
        return

    def get_deform_local(self, node_set: np.ndarray) -> list:
        """
        :param node_set:
            节点信息

        :return:
            节点的坐标信息
        """
        deform_local = []
        for node in node_set:
            u = self.result_dict['deform']['Ux'][node]
            v = self.result_dict['deform']['Uy'][node]
            deform_local.append(u)
            deform_local.append(v)
        return deform_local

    def post_process(self):
        """
        根据计算的结果获取各节点的应力信息
        """
        self.strain_list = []
        self.stress_list = []
        for element_index, node_set in enumerate(self.elements):
            element = self.instant_element(element_index, node_set)
            deform_local = self.get_deform_local(node_set)
            strain = element.get_strain(deform_local)
            self.strain_list.append(strain)
            stress = element.get_strain(deform_local)
            self.stress_list.append(stress)
        self.strain_list = np.array(self.strain_list).T
        self.stress_list = np.array(self.stress_list).T
        return

    def update_show_dict(self):
        """
        根据计算的结果和计算得到的应力信息，
        更新结果字典
        """
        self.result_dict['deform']['Ux'] = self.deform.reshape(len(self.nodes), 2).T[0]
        self.result_dict['deform']['Uy'] = self.deform.reshape(len(self.nodes), 2).T[1]
        self.result_dict['force']['Fx'] = self.force.reshape(len(self.nodes), 2).T[0]
        self.result_dict['force']['Fy'] = self.force.reshape(len(self.nodes), 2).T[1]
        self.post_process()
        for i, name in enumerate(('e11', 'e22', 'e12')):
            self.result_dict['strain'][name] = self.strain_list[i]
        for i, name in enumerate(('S11', 'S22', 'S12')):
            self.result_dict['stress'][name] = self.stress_list[i]
        return

    def solve(self):
        """
        首先获取待求解的 FEM 方程，
        随后使用 numpy 进行求解，
        最后更新相应的变量
        """
        self.get_fem()
        self.solve_fem()
        self.update_global_variables()
        self.update_show_dict()
        return

    def plot_mesh(self, add_text: bool):
        """
        绘制三角形网格

        :param add_text:
            是否添加节点编号信息
        """
        _x = self.result_dict['position']['x']
        _y = self.result_dict['position']['y']
        mesh_fig = tri.Triangulation(_x, _y, self.elements)
        plt.figure()
        plt.gca().set_aspect('equal')
        plt.triplot(mesh_fig, 'b.-', lw=1)
        nodes = self.mesh.points
        if add_text:
            for i, _ in enumerate(nodes):
                plt.annotate(str(i), (_x[i], _y[i]), color="red")
        plt.title('三角形网格')
        plt.show()
        return

    def plot_result(self, name, component, shading):
        """
        汇总 FEM 的计算结果可视化图
        """
        _x = self.result_dict['position']['x']
        _y = self.result_dict['position']['y']
        _z = self.result_dict[name][component]
        fig = tri.Triangulation(_x, _y, self.elements)
        plt.figure()
        plt.gca().set_aspect('equal')
        plt.tripcolor(fig, _z, shading=shading)
        plt.triplot(fig, lw=1)
        plt.colorbar()
        plt.title(f'FEM 计算结果: {name}')
        plt.show()
        return


if __name__ == '__main__':
    import pygmsh
    import matplotlib as mpl

    # 将字体设置为 等线，
    # 以便于显示中文文本，
    mpl.rc("font", family='DengXian')
    mpl.rcParams['axes.unicode_minus'] = False

    # 使用 pygmsh 生成三角形网格，不参与计算
    with pygmsh.geo.Geometry() as geom:
        geom.add_polygon(
            [
                [0.0, 0.0],
                [2.0, 0.0],
                [2.0, 1.0],
                [0.0, 1.0],
            ],
            mesh_size=0.1,
        )
        mesh = geom.generate_mesh()

    # 设置材料
    material_test = Material(e=10.0, nv=0.5)
    print(material_test.D)
    Fea.set_material(mesh, material_test)

    fea = Fea(mesh)

    # 设置边界约束以及外加载荷
    # 其中约束为单边固定
    # 载荷是单边中点处施加平行力
    x_fix_test = {}
    y_fix_test = {}
    f_given_test = {}
    for index, position in enumerate(mesh.points):
        x = position[0].tolist()
        y = position[1].tolist()
        if x < 1e-6:  # 设置为 0 时因计算精度问题会出现错误
            x_fix_test.update({index: 0.0})
            y_fix_test.update({index: 0.0})
        if 1.0 - x < 1e-6 and abs(y - 0.5) < 1e-6:
            f_given_test.update({index: (1.0, 0.0)})  # 外加载荷的大小为 Fx=1 Fy=0

    # 提交边界信息并计算
    fea.set_conditions(x_fix_test, y_fix_test, f_given_test)
    fea.solve()

    # 汇总网格可视化图
    fea.plot_mesh(False)
    fea.plot_mesh(True)

    # 绘制形变和应力可视化图
    fea.plot_result('deform', 'Ux', 'gouraud')
    fea.plot_result('strain', 'e11', 'flat')
