import numpy as np

from wzc.utilities import tools as utils
from wzc.fem.fem2d import FEM


class FEMSetup(FEM):

    def __init__(
            self,
            young_modulus,
            poisson_ratio,
            **kwargs
    ):
        super(FEMSetup, self).__init__(**kwargs)
        self.young_modulus = young_modulus
        self.poisson_ratio = poisson_ratio

        self.lamb = (self.young_modulus * self.poisson_ratio) / ((1 + self.poisson_ratio) * (1 - 2 * self.poisson_ratio))
        self.mu = self.young_modulus / (2 * (1 + self.poisson_ratio))
        return

    @staticmethod
    def stima3aux_symetric(coords: np.ndarray):
        """
        方程中的对称项部分的计算
        """
        t = utils.area_of_triangle(coords)

        x_diff = np.roll(coords[:, 0], 1) - np.roll(coords[:, 0], 2)
        x_diff_mat = np.expand_dims(x_diff, 1) @ np.expand_dims(x_diff, 0)

        y_diff = np.roll(coords[:, 1], -1) - np.roll(coords[:, 1], -2)
        y_diff_mat = np.expand_dims(y_diff, 1) @ np.expand_dims(y_diff, 0)

        return t, x_diff_mat, y_diff_mat

    @staticmethod
    def stima3aux_asymetric(coords: np.ndarray):
        """
        方程中的非对称项部分的计算
        """
        t = utils.area_of_triangle(coords)

        x_diff = np.roll(coords[:, 0], 1) - np.roll(coords[:, 0], 2)
        y_diff = np.roll(coords[:, 1], -1) - np.roll(coords[:, 1], -2)
        xy_diff_mat = np.expand_dims(x_diff, 1) @ np.expand_dims(y_diff, 0)

        return t, xy_diff_mat

    def stima3xx(self, coords: np.ndarray) -> np.ndarray:
        t, x_diff_mat, y_diff_mat = self.stima3aux_symetric(coords)
        return ((self.mu + self.lamb) * y_diff_mat + 0.5 * self.mu * x_diff_mat) / (4 * t)

    def stima3yy(self, coords: np.ndarray) -> np.ndarray:
        t, x_diff_mat, y_diff_mat = self.stima3aux_symetric(coords)
        return ((self.mu + self.lamb) * x_diff_mat + 0.5 * self.mu * y_diff_mat) / (4 * t)

    def stima3xy(self, coords: np.ndarray) -> np.ndarray:
        t, xy_diff_mat = self.stima3aux_asymetric(coords)
        return (self.lamb * xy_diff_mat + 0.5 * self.mu * xy_diff_mat.T) / (4 * t)

    def stima3yx(self, coords: np.ndarray) -> np.ndarray:
        t, xy_diff_mat = self.stima3aux_asymetric(coords)
        return (self.lamb * xy_diff_mat.T + 0.5 * self.mu * xy_diff_mat) / (4 * t)

    def solve(self):
        """
        根据泊松方程对线弹性方程求解
        """
        nodes_num = self.mesh.nodes_num
        base_func_num = 2 * nodes_num
        a = np.zeros(shape=(base_func_num, base_func_num))
        b = np.zeros(shape=(base_func_num,))

        for func, beg_x, beg_y in [
            (self.stima3xx, 0, 0),
            (self.stima3xy, nodes_num, 0),
            (self.stima3yx, 0, nodes_num),
            (self.stima3yy, nodes_num, nodes_num)
        ]:
            for nodes in self.mesh.nodes_of_elem:
                local = func(self.mesh.coordinates2D[nodes])
                for y, col in enumerate(nodes):
                    for x, row in enumerate(nodes):
                        a[row + beg_x, col + beg_y] += local[x, y]

        for nodes in self.mesh.nodes_of_elem:
            rhs_x, rhs_y = self.rhs_val(self.mesh.coordinates2D[nodes])
            b[nodes] += rhs_x
            b[nodes + nodes_num] += rhs_y

        for vert_idxs in self.mesh.neumann_boundaries:
            neu_x, neu_y = self.neumann_val(self.mesh.coordinates2D[vert_idxs])
            b[vert_idxs] += neu_x
            b[vert_idxs + nodes_num] += neu_y

        u = np.zeros(shape=(base_func_num, 1))
        for vert_idxs in self.mesh.dirichlet_boundaries:
            values = self.dirichlet_val(self.mesh.coordinates2D[vert_idxs])
            dir_x, dir_y = values[:, 0], values[:, 1]
            u[vert_idxs] = np.expand_dims(dir_x, 1)
            u[vert_idxs + nodes_num] = np.expand_dims(dir_y, 1)
        b -= (a @ u).T[0]

        free_nodes = [v for v in range(nodes_num) if v not in self.mesh.dirichlet_boundaries]
        free_nodes += [v + nodes_num for v in free_nodes]

        u_free = np.linalg.solve(a[free_nodes][:, free_nodes], b[free_nodes])
        w = np.squeeze(u.copy(), 1)
        w[free_nodes] = u_free
        return w
