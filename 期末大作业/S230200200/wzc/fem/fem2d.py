from typing import Optional, Callable, Union
import numpy as np

from wzc.utilities import tools as utils
from wzc.mesh import Mesh


class FEM:
    """
    FEM 父类
    """
    def __init__(
            self,
            mesh: Mesh,
            rhs_func: Callable[[np.ndarray], Union[float, np.ndarray]],
            dirichlet_func: Optional[Callable[[np.ndarray], Union[float, np.ndarray]]] = None,
            neumann_func: Optional[Callable[[np.ndarray], Union[float, np.ndarray]]] = None
    ):
        self.mesh = mesh
        self.rhs_func = rhs_func
        self.dirichlet_func = dirichlet_func
        self.neumann_func = neumann_func
        return

    def rhs_val(self, vertices: np.ndarray) -> float:
        t = utils.area_of_triangle(vertices)
        center = utils.center_of_mass(vertices)
        return t / 3 * self.rhs_func(center)

    def neumann_val(self, vertices: np.array) -> float:
        if vertices.shape[0] == 1:
            return self.neumann_func(vertices)
        center = utils.center_of_mass(vertices)
        length = np.sqrt((vertices[1, 0] - vertices[0, 0]) ** 2 + (vertices[1, 1] - vertices[0, 1]) ** 2)

        return self.neumann_func(center) * length / 2

    def dirichlet_val(self, vertices: np.array) -> np.ndarray:
        values = np.array([self.dirichlet_func(v) for v in vertices])
        return values

    def solve(self):
        raise NotImplementedError
