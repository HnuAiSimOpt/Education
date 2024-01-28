from enum import Enum
from typing import List, Dict

import numpy as np
import meshio


class Mesh:
    """
    参考 GMSH 官方的读取示例

    读取标准 GMSH 划分后的网格
    """

    # 类常量
    class ElementType(Enum):
        VERTEX = 'vertex'
        LINE = 'line'
        TRIANGULAR = 'triangle'

    class BoundaryConditionType(Enum):
        DIRICHLET = 'dirichlet'
        NEUMANN = 'neumann'

    def __init__(self, mesh_filename: str):

        msh = meshio.read(mesh_filename)

        self.nodes_of_elem = msh.cells_dict[self.ElementType.TRIANGULAR.value]
        self.coordinates = msh.points
        self.coordinates2D = self.coordinates[:, 0:2]
        self.nodes_num = self.coordinates.shape[0]
        self.elems_num = len(self.nodes_of_elem)

        self.physical_groups_mapping = self.extract_physical_groups(msh)
        self.dirichlet_boundaries = np.array([])
        self.neumann_boundaries = np.array([])
        return

    def extract_physical_groups(self, msh) -> Dict[str, np.ndarray]:
        # 读取预留的材料位
        condition_of_elem = self.prepare_physical_groups_mapping(msh)
        physical_groups_mapping = {}

        physical_groups = msh.cell_data['gmsh:physical']
        group_idx = 0
        inner_idx = 0
        for elem_type, elem_nodes in msh.cells_dict.items():

            for node in elem_nodes:
                key = physical_groups[group_idx][inner_idx]
                condition = condition_of_elem[key]
                if condition not in physical_groups_mapping:
                    physical_groups_mapping[condition] = []
                physical_groups_mapping[condition].append(node)

                inner_idx += 1
                if inner_idx >= physical_groups[group_idx].size:
                    inner_idx = 0
                    group_idx += 1

        return {k: np.array(v) for k, v in physical_groups_mapping.items()}

    @staticmethod
    def prepare_physical_groups_mapping(msh) -> list:
        condition_of_elem = ['undefined'] * (max([v[0] for v in msh.field_data.values()]) + 1)
        for group_name, val in msh.field_data.items():
            condition_of_elem[val[0]] = group_name
        return condition_of_elem

    def set_boundary_condition(self, boundary_type: BoundaryConditionType, groups: List[str]):
        # 设置边界和施加的约束
        boundaries = np.concatenate(tuple(self.physical_groups_mapping[name] for name in groups))
        if boundary_type == Mesh.BoundaryConditionType.DIRICHLET:
            self.dirichlet_boundaries = boundaries
        elif boundary_type == Mesh.BoundaryConditionType.NEUMANN:
            self.neumann_boundaries = boundaries
        return
