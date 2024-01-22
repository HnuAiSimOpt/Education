import numpy as np
import scipy as s
from PIL import Image

"""
   本代码为求指定三角形单元和各节点在给定载荷F下的应力和位移
"""

image = Image.open('question.jpg')
image.show()

"""
   定义全局刚度矩阵和载荷矩阵
"""

global stiffness_matrix
global F

"""
   需要用到的函数
"""

"""
   生成荷载向量F
"""
def load_F(P_inp, num_nodes):

    F = np.zeros((2 * num_nodes, 1))
    for node_load in P_inp:
        node_index = int(node_load[0]) - 1
        F[2 * node_index] = node_load[1]    # X方向荷载
        F[2 * node_index + 1] = node_load[2]  # Y方向荷载
    return F

"""
   处理边界条件
"""

def apply_boundary_conditions(stiffness_matrix, F, boundary):

    for k, (node, ux_bc, uy_bc) in enumerate(boundary.astype(int)):
        i = 2 * (node - 1)
        if ux_bc == 1:
            F[i] = 0
            stiffness_matrix[i, :] = 0
            stiffness_matrix[:, i] = 0
            stiffness_matrix[i, i] = 1
        if uy_bc == 1:
            F[i + 1] = 0
            stiffness_matrix[i + 1, :] = 0
            stiffness_matrix[:, i + 1] = 0
            stiffness_matrix[i + 1, i + 1] = 1
    return stiffness_matrix, F

"""
   计算节点位移
"""

def get_nodal_displacements(U, num_nodes):

    Uout = np.zeros((num_nodes, 3))
    Uout[:, 0] = np.arange(1, num_nodes + 1)
    Uout[:, 1] = U[0::2]
    Uout[:, 2] = U[1::2]

    return Uout

"""
   计算单元应力
"""

def solve_stress(U,nodes, elements, num_elements, elastic_modulus, poisson_ratio):

    stress_matrix = np.zeros((num_elements, 4))

    for element_index in range(num_elements):
        displacement_vector = np.zeros((6, 1))

        element_nodes = elements[element_index, :3] - 1
        displacement_vector[::2, 0] = U[element_nodes * 2]
        displacement_vector[1::2, 0] = U[element_nodes * 2 + 1]

        node_coords = nodes[element_nodes, :2].flatten()
        x_coords, y_coords = node_coords[::2], node_coords[1::2]
        x_extended, y_extended = np.append(x_coords, x_coords[:2]), np.append(y_coords, y_coords[:2])

        a = x_extended[1:4] * y_extended[2:5] - x_extended[2:5] * y_extended[1:4]
        b = y_extended[1:4] - y_extended[2:5]
        c = -x_extended[1:4] + x_extended[2:5]

        matrix_a = np.vstack(([1, 1, 1], x_coords, y_coords)).T
        area = 0.5 * np.linalg.det(matrix_a)

        stress_element_matrix = np.zeros((3, 6))

        elastic_term = elastic_modulus / (2 * (1 - poisson_ratio ** 2) * area)
        for i in range(3):
            stress_element_matrix[:, 2 * i:2 * (i + 1)] = \
                elastic_term * np.array([[b[i], poisson_ratio * c[i]],
                                         [poisson_ratio * b[i], c[i]],
                                         [(1 - poisson_ratio) / 2 * c[i], (1 - poisson_ratio) / 2 * b[i]]])

        element_stress = stress_element_matrix @ displacement_vector #计算应力
        stress_matrix[element_index, :] = [element_index + 1, *element_stress.flatten()]

    return stress_matrix

"""
   生成刚度矩阵
"""

def stiff_matrix( E, poisson_ratio, T, num_nodes,num_elements, nodes, elements):

   Ns = 2 * num_nodes
   stiffness_matrix = np.zeros((Ns,Ns))

   for n in range(1,num_elements+1):
       x = np.array([0] * 6);y = np.array([0] * 6)

       for i in range(3):
           x[i] = nodes[elements[n - 1, i] - 1, 0]
           y[i] = nodes[elements[n - 1, i] - 1, 1]

       x[3:6] = x[0:3]
       y[3:6] = y[0:3]

       MA = np.array([[1,x[0],y[0]],[1,x[1],y[1]],[1,x[2],y[2]]])
       A = 0.5*np.linalg.det(MA)

       a = np.array([0] * 3);
       b = np.array([0] * 3);
       c = np.array([0] * 3)

       for i in range(3):
           a[i] = x[i + 1] * y[i + 2] - x[i + 2] * y[i+1]
           b[i] = y[i+1] - y[i + 2]
           c[i] = -x[i+1] + x[i + 2]

       EK = np.zeros((6,6))

       for i in range(3):
           for j in range(3):

               k1 = b[i]*b[j]+(1-poisson_ratio)/2*c[i]*c[j]
               k2 = poisson_ratio*c[i]*b[j]+(1-poisson_ratio)/2*b[i]*c[j]
               k3 = poisson_ratio*b[i]*c[j]+(1-poisson_ratio)/2*c[i]*b[j]
               k4 = c[i]*c[j]+(1-poisson_ratio)/2*b[i]*b[j]
               EK[2*i:2*i+2,2*j:2*j+2] = E*T/(4*(1-poisson_ratio**2)*A)*np.array([[k1,k3],[k2,k4]])

       k_temp = np.zeros((Ns,Ns))
       temp = np.zeros((6,Ns))

       for i in range(3):
           temp[2*i,2*elements[n-1,i]-2] = 1
           temp[2*i+1,2*elements[n -1,i]-1] = 1

       temp_T = np.transpose(temp)
       k_temp=np.dot(np.dot(temp_T,EK),temp)
       stiffness_matrix = stiffness_matrix + k_temp

   return stiffness_matrix

"""
   主函数
"""

"""
   输入相关参数
"""

E = 20
poisson_ratio = 0
T = 1
num_nodes = 4
num_elements = 2
nodes = np.array([[0,0], [0,1], [2,1], [1,0]])
elements = np.array([[1,2,3], [3,4,1]])
P_inp= np.array([[4,0,-5]])
boundary=np.array([[1,1,1],[2,1,1]])

"""
  求节点位移
"""

stiffness_matrix=stiff_matrix(E, poisson_ratio, T, num_nodes, num_elements, nodes, elements)
F=load_F(P_inp,num_nodes)

apply_boundary_conditions(stiffness_matrix,F,boundary)
U=s.sparse.linalg.spsolve(stiffness_matrix,F)

"""
   求应力、应变
"""

Uout=get_nodal_displacements(U,num_nodes)
S=solve_stress(U,nodes,elements,num_elements,E,poisson_ratio)

"""
   输出结果
"""
print('各节点位移:'
      '(第一列为节点编号，第二列第三列分别为节点X、Y方向位移)')
print(Uout)
print('各单元应力:'
      '(第一列为单元编号，第二列、第三列和第四列分别为 σx, σy, τxy)')
print(S)