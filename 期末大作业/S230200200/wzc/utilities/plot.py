import numpy as np
from matplotlib import pyplot as plt, tri

from wzc.mesh import Mesh


def plot_results(mesh: Mesh, magnitudes: np.ndarray, file_name: str = 'results.png'):
    """
    绘制结果曲线
    """
    triangulation = tri.Triangulation(
        x=mesh.coordinates2D[:, 0],
        y=mesh.coordinates2D[:, 1],
        triangles=mesh.nodes_of_elem
    )
    plt.tricontourf(triangulation, magnitudes)
    plt.colorbar()
    plt.savefig(file_name)
    plt.close()
    return


def plot_displacements(mesh: Mesh, displacements: np.ndarray, file_name: str = 'displacements.png'):
    """
    绘制形变
    """
    before = tri.Triangulation(
        x=mesh.coordinates2D[:, 0],
        y=mesh.coordinates2D[:, 1],
        triangles=mesh.nodes_of_elem
    )
    plt.triplot(before, color='#1f77b4')
    after = tri.Triangulation(
        x=mesh.coordinates2D[:, 0] + displacements[:, 0],
        y=mesh.coordinates2D[:, 1] + displacements[:, 1],
        triangles=mesh.nodes_of_elem
    )
    plt.triplot(after, color='#ff7f0e')
    plt.grid()
    plt.savefig(file_name)
    plt.close()
    return


def plot_stress(mesh: Mesh, stress: np.ndarray, file_name: str = 'stress.png'):
    """
    绘制应力
    """
    triangulation = tri.Triangulation(
        x=mesh.coordinates2D[:, 0],
        y=mesh.coordinates2D[:, 1],
        triangles=mesh.nodes_of_elem
    )
    plt.tripcolor(triangulation, stress)
    plt.colorbar()
    plt.savefig(file_name)
    plt.close()
    return
