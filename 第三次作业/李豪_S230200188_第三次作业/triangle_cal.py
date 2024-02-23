# 姓名：李豪  学号：S230200188
# 问题假定了一块矩形薄板，长0.5米，宽0.025米，高0.25米；
# 薄板左边固定，右边承受了均匀分布的拉力3000kN/m^2,即3Mpa的应力
# 四节点坐标逆时针排列
# (x1,y1)=(0,0)       (x2,y2)=(0.5,0)  
# (x3,y3)=(0.5,0.25)  (x4,y4)=(0,0.25)
# 运行后结果: 四节点的受力为变量F，四节点的位移为变量U，
# 单元1和2的应力为sigma1和sigma2,应变为epsilon1和epsilon2
import numpy as np


def Linear_Triangle_Element_Stiffness(E,NU,t,xi,yi,xj,yj,xm,ym,p):
    #弹性模量为 E、泊松比为 NU、厚度为t 
    #p=1用于面应力情况 p=2表明函数用于平面应变情况。
    #该函数返回6X6的单元刚度矩阵。
    A =(xi*(yj-ym) + xj*(ym-yi) + xm*(yi-yj))/2
    betai = yj-ym
    betaj = ym-yi
    betam = yi-yj
    gammai = xm-xj
    gammaj = xi-xm
    gammam = xj-xi
    B = np.array([[betai, 0, betaj, 0, betam, 0],
                   [0, gammai, 0, gammaj, 0, gammam],
                   [gammai, betai, gammaj, betaj, gammam, betam]])/(2*A)
    if p==1:
        D=(E/(1-NU*NU))*np.array([[1, NU, 0],
                                [NU, 1, 0],
                                [0, 0, (1-NU)/2]])
    elif p == 2:
        D=(E/(1+NU)/(1-2*NU) )*np.array([[1-NU, NU, 0],
                                        [NU, 1-NU, 0],
                                         [0, 0, (1-2*NU)/2]])
    # 计算 B 的转置
    B_transpose = np.transpose(B)
    # 执行 B' * D * B
    return t*A*np.dot(np.dot(B_transpose, D), B)


def Linear_Triangle_Assemble(K,k,i,j,m):
    #该函数将连接节点 i和节点j的线性三角形元的
    #单元刚度矩阵k集成到整体刚度矩阵K。
    #每集成一个单元，该函数都将返 2nX2n 的整体刚度矩阵K。

    K[2*i-2, 2*i-2] = K[2*i-2, 2*i-2] + k[0, 0]
    K[2*i-2, 2*i-1] = K[2*i-2, 2*i-1] + k[0, 1]
    K[2*i-2, 2*j-2] = K[2*i-2, 2*j-2] + k[0, 2]
    K[2*i-2, 2*j-1] = K[2*i-2, 2*j-1] + k[0, 3]
    K[2*i-2, 2*m-2] = K[2*i-2, 2*m-2] + k[0, 4]
    K[2*i-2, 2*m-1] = K[2*i-2, 2*m-1] + k[0, 5]

    K[2*i-1, 2*i-2] = K[2*i-1, 2*i-2] + k[0, 0]
    K[2*i-1, 2*i-1] = K[2*i-1, 2*i-1] + k[0, 1]
    K[2*i-1, 2*j-2] = K[2*i-1, 2*j-2] + k[0, 2]
    K[2*i-1, 2*j-1] = K[2*i-1, 2*j-1] + k[0, 3]
    K[2*i-1, 2*m-2] = K[2*i-1, 2*m-2] + k[0, 4]
    K[2*i-1, 2*m-1] = K[2*i-1, 2*m-1] + k[0, 5]

    K[2*j-2, 2*i-2] = K[2*j-2, 2*i-2] + k[2, 0]
    K[2*j-2, 2*i-1] = K[2*j-2, 2*i-1] + k[2, 1]
    K[2*j-2, 2*j-2] = K[2*j-2, 2*j-2] + k[2, 2]
    K[2*j-2, 2*j-1] = K[2*j-2, 2*j-1] + k[2, 3]
    K[2*j-2, 2*m-2] = K[2*j-2, 2*m-2] + k[2, 4]
    K[2*j-2, 2*m-1] = K[2*j-2, 2*m-1] + k[2, 5]

    K[2*j-1, 2*i-2] = K[2*j-1, 2*i-2] + k[3, 0]
    K[2*j-1, 2*i-1] = K[2*j-1, 2*i-1] + k[3, 1]
    K[2*j-1, 2*j-2] = K[2*j-1, 2*j-2] + k[3, 2]
    K[2*j-1, 2*j-1] = K[2*j-1, 2*j-1] + k[3, 3]
    K[2*j-1, 2*m-2] = K[2*j-1, 2*m-2] + k[3, 4]
    K[2*j-1, 2*m-1] = K[2*j-1, 2*m-1] + k[3, 5]

    K[2*m-2, 2*i-2] = K[2*m-2, 2*i-2] + k[4, 0]
    K[2*m-2, 2*i-1] = K[2*m-2, 2*i-1] + k[4, 1]
    K[2*m-2, 2*j-2] = K[2*m-2, 2*j-2] + k[4, 2]
    K[2*m-2, 2*j-1] = K[2*m-2, 2*j-1] + k[4, 3]
    K[2*m-2, 2*m-2] = K[2*m-2, 2*m-2] + k[4, 4]
    K[2*m-2, 2*m-1] = K[2*m-2, 2*m-1] + k[4, 5]

    K[2*m-1, 2*i-2] = K[2*m-1, 2*i-2] + k[5, 0]
    K[2*m-1, 2*i-1] = K[2*m-1, 2*i-1] + k[5, 1]
    K[2*m-1, 2*j-2] = K[2*m-1, 2*j-2] + k[5, 2]
    K[2*m-1, 2*j-1] = K[2*m-1, 2*j-1] + k[5, 3]
    K[2*m-1, 2*m-2] = K[2*m-1, 2*m-2] + k[5, 4]
    K[2*m-1, 2*m-1] = K[2*m-1, 2*m-1] + k[5, 5]
    return K


def Linear_Triangle_Element_Stresses(E, NU, xi, yi, xj, yj, xm, ym, p, u):    
    #该函数计算在LinearTriangleElementStresses(E,NU,t,x y xp yp xm,ym p,u)
    #弹性模量为 E、泊松比为 NU、以及单元位移矢量为u的单元应力。
    #p=1 用于面应力情况。p=2用于平面应变情况。
    #该函数返回单元应力矢量。
    A=(xi*(yj-ym)+ xj*(ym-yi) + xm*(yi-yj))/2
    betai = yj-ym
    betaj = ym-yi
    betam = yi-yj
    gammai = xm-xj
    gammaj = xi-xm
    gammam = xj-xi
    B =np.array([[betai, 0, betaj, 0, betam, 0],
                   [0, gammai, 0, gammaj, 0, gammam],
                   [gammai, betai, gammaj, betaj, gammam, betam]])/(2*A)
    if p==1:
        D=(E/(1-NU*NU))*np.array([[1, NU, 0],
                                [NU, 1, 0],
                                [0, 0, (1-NU)/2]])
    elif p == 2:
        D=(E/(1+NU)/(1-2*NU) )*np.array([[1-NU, NU, 0],
                                        [NU, 1-NU, 0],
                                         [0, 0, (1-2*NU)/2]])
    return D @ B @ u
        

def Linear_Triangle_Element_PStresses(sigma):
    #该函数根据单元应力矢 sigma 计算 单元主应力。
    #LinearTriangleElementStresses(sigma)该函数
    #返回3X1的矢量,其形式为[sigmal,sigma2,theta]。其中 sigmal 和 sigma2为单元的主应力,
    #theta 为主应力方向角。 
    R =(sigma[0] + sigma[1])/2
    Q =((sigma[0] - sigma[1])/2)**2 + sigma[2]*sigma[2]
    M = 2 * np.divide(sigma[2], (sigma[0] - sigma[1]))
    s1= R + np.sqrt(Q)
    s2= R - np.sqrt(Q)
    theta = (np.arctan(M)/2)*180/np.pi
    return(s1,s2,theta)


E=210e6 #弹性模量 单位Pa
NU=0.3  #泊松比
t=0.025 #薄板厚度 单位m

# 计算求解
# 计算第一个单元刚度矩阵（1,3,4）节点
k1=Linear_Triangle_Element_Stiffness(E,NU,t,0,0,0.5,0.25,0,0.25,1)
#计算第二个单元刚度矩阵（1,2,3）节点
k2=Linear_Triangle_Element_Stiffness(E,NU,t,0,0,0.5,0,0.5, 0.25,1)
#生成总体刚度矩阵
K=np.zeros((8,8))
K=Linear_Triangle_Assemble(K,k1,1,3,4)
K=Linear_Triangle_Assemble(K,k2,1,2,3)
#求出节点位移和受力
k=K[2:6, 2:6]
f=np.array([[9.375],[0],[9.375],[0]])
u=np.linalg.solve(k, f)   
U=np.zeros((4,1))
U=np.insert(U, 2, u, axis=0)
F=np.dot(K,U)

#求解单元应力应变
u1=np.array([U[0],U[1],U[4],U[5],U[6],U[7]])
u2=np.array([U[0],U[1],U[2],U[3],U[4],U[5]])
np.set_printoptions(precision=3)
print('-'*40)
print('四个节点受力F=\n',F.reshape(4,2)*1000,'\nN',)
print('-'*40)
print('单元1（1,3,4）的节点位移为u1=\n',u1.reshape(3,2)*1e5,'\ncm','\n','-'*40,'\n'
      '单元2（1,2,3）的节点位移为u2=\n',u2.reshape(3,2)*1e5,'\ncm')
sigma1=Linear_Triangle_Element_Stresses(E,NU,0,0,0.5,0.25,0,0.25,1,u1)
epsilon1=Linear_Triangle_Element_Stresses(E,NU,0,0,0.5,0.25,0,0.25,2,u1)
sigma2=Linear_Triangle_Element_Stresses(E,NU,0,0,0.5,0,0.5,0.25,1,u2)
epsilon2=Linear_Triangle_Element_Stresses(E,NU,0,0,0.5,0,0.5,0.25,2,u2)
s1=Linear_Triangle_Element_PStresses(sigma1)
s2=Linear_Triangle_Element_PStresses(sigma2)
print('-'*40)
print('单元1（1,3,4）的应力为sigma1=\n',sigma1.reshape(1,3)*1e-3,'\nMpa','\n','-'*40,'\n'
      '单元2（1,2,3）的应力为sigma2=\n',sigma2.reshape(1,3)*1e-3,'\nMpa')
print('-'*40)
 
  
    