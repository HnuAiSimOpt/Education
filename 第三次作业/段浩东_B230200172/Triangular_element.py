"""
   说明：
   作者：段浩东
   学号：B230200172
   程序目的：通过三角形单元求解一个5单元7节点的算例，得到每一个节点的位移。具体内容见附带的图片文件。
"""

import numpy as np
import scipy as s
"""
   ###子函数部分###
"""
def create_K( EL, PO, TK, NODE_SUM,ELE_SUM, NODE, ELE):  # 生成整体刚度矩阵
   """
   生成单元刚度矩阵并添加到整体刚度矩阵"ASK"
   PO：泊松比  NODE_SUM：节点数
   EL：弹性模量 TK：厚度 NODE：节点坐标矩阵   ELE：单元节点号矩阵
   返回值：整体刚度矩阵
   """
   ASK = np.zeros((2 * NODE_SUM, 2 * NODE_SUM))  # 整体刚度矩阵初始化（ASK为局部变量）
   for n in range(1,ELE_SUM+1):
       ##1.整理当前单元节点坐标##
       x = np.array([0] * 6);y = np.array([0] * 6)   # 当前单元坐标数组
       for i in range(3):  # 将当前单元的坐标存入数组，顺序为单元节点号矩阵内的顺序，编号从0开始(x0,x1,x2)
           x[i] = NODE[ELE[n - 1, i] - 1, 0]
           y[i] = NODE[ELE[n - 1, i] - 1, 1]
       x[3:6] = x[0:3]
       y[3:6] = y[0:3]     # 便于后续下标轮换调用
       ##2.计算三角形单元面积A
       MA=np.array([[1,x[0],y[0]],[1,x[1],y[1]],[1,x[2],y[2]]])
       A=0.5*np.linalg.det(MA)
       ##3.计算插值函数系数
       a = np.array([0] * 3);
       b = np.array([0] * 3);
       c = np.array([0] * 3)  # a,b,c为插值函数系数，初始化
       for i in range(3):
           a[i] = x[i + 1] * y[i + 2] - x[i + 2] * y[i+1]  # ai=xjym-xmyj,下标i，j，m循环(a0,a1,a2)
           b[i] = y[i+1] - y[i + 2]  # bi=yi-ym
           c[i] = -x[i+1] + x[i + 2]  # ci=-xi+xm
       ##4.生成单元刚度矩阵##
       EK=np.zeros((6,6))
       for i in range(3):
           for j in range(3):
               k1=b[i]*b[j]+(1-PO)/2*c[i]*c[j]
               k2=PO*c[i]*b[j]+(1-PO)/2*b[i]*c[j]
               k3=PO*b[i]*c[j]+(1-PO)/2*c[i]*b[j]
               k4=c[i]*c[j]+(1-PO)/2*b[i]*b[j]
               EK[2*i:2*i+2,2*j:2*j+2]=EL*TK/(4*(1-PO**2)*A)*np.array([[k1,k3],[k2,k4]])
       ##4.将单元刚度矩阵组装为整体刚度矩阵
       ask=np.zeros((2*NODE_SUM,2*NODE_SUM)) #临时总刚初始化
       G=np.zeros((6,2*NODE_SUM))            #定义转换矩阵G
       for i in range(3):
           G[2*i,2*ELE[n-1,i]-2]=1
           G[2*i+1,2*ELE[n -1,i]-1]=1
       GT=np.transpose(G)
       ask=np.dot(np.dot(GT,EK),G)
       ASK=ASK+ask
   return ASK


def create_P(P_inp,NODE_SUM):#由输入的节点荷载信息生成荷载向量P
   """
   由输入的节点荷载信息生成荷载向量P
   P_inp:节点荷载矩阵   NODE_SUM：节点数
   返回值：荷载向量
   """
   P=np.zeros((2*NODE_SUM,1))  #初始化荷载向量P
   k=1
   P2 = P_inp.astype(int)  # 确保循环中的i为整型
   for i in P2[:,0]:
       P[2*i-2]=P_inp[k-1,1]
       P[2*i-1]=P_inp[k-1,2]
       k=k+1
   return (P)

def do_BC(BC):#处理边界条件，置1法对角元素
   """
   修改主函数中总体刚度矩阵ASK，和荷载向量P
   """
   global ASK  #申明全局变量，以便修改
   global P
   k=1
   for i in BC[:,0]:
       if BC[k-1,1]==1:
           P[2*i-2]=0
       if BC[k-1,2]==1:
           P[2*i-1]=0
       k=k+1
   k=1
   for i in BC[:, 0]:
       if BC[k-1, 1]==1:
           ASK[2*i-2,:]=0
           ASK[:,2*i-2]=0
           ASK[2*i-2,2*i-2]=1
       if BC[k-1,2]==1:
           ASK[2*i-1,:] = 0
           ASK[:,2*i-1] = 0
           ASK[2*i-1,2*i-1]=1
       k=k+1
def solve_Uout(U,NODE_SUM):
   """
   #根据位移向量得到节点位移信息矩阵
   返回值：结构应变信息矩阵Uout
   第一列：单元号 第二列~第三列：Ux,Uy
   """
   Uout=np.zeros((NODE_SUM,3))
   for i in range(NODE_SUM):
       Uout[i,0]=i+1
       Uout[i,1]=U[2*i]
       Uout[i, 2]=U[2*i+1]
   return Uout

def solve_B(U,NODE,ELE,ELE_SUM):
   """
   #根据位移向量U求解结构应变
   返回值：结构应变信息矩阵B
   第一列：单元号 第二列~第四列：εx,εy,γxy
   """
   B = np.zeros((ELE_SUM, 4))
   for n in range(1, ELE_SUM + 1):
       ##1.建立当前单元的位移向量
       ue = np.zeros((6, 1))  # 初始化单元位移向量ue
       for i in range(3):
           ue[2 * i, 0] = U[2 * ELE[n - 1, i] - 2]
           ue[2 * i + 1, 0] = U[2 * ELE[n - 1, i] - 1]
       ##2.整理当前单元节点坐标##
       x = np.array([0] * 6)
       y = np.array([0] * 6)  # 当前单元坐标数组
       for i in range(3):     # 将当前单元的坐标存入数组x,y，顺序为单元节点号矩阵内的顺序，编号从0开始(x0,x1,x2)
           x[i] = NODE[ELE[n - 1, i] - 1, 0]
           y[i] = NODE[ELE[n - 1, i] - 1, 1]
       x[3:6] = x[0:3]
       y[3:6] = y[0:3]        # 便于后续下标轮换调用
       ##3.计算插值函数系数
       a = np.array([0] * 3);
       b = np.array([0] * 3);
       c = np.array([0] * 3)  # a,b,c为插值函数系数，初始化
       for i in range(3):
           a[i] = x[i + 1] * y[i + 2] - x[i + 2] * y[i+1]  # ai=xjym-xmyj,下标i，j，m循环(a0,a1,a2)
           b[i] = y[i + 1] - y[i + 2]  # bi=yi-ym
           c[i] = -x[i + 1] + x[i + 2]  # ci=-xi+xm
       ##4.计算三角形单元面积A
       MA = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])
       A = 0.5 * np.linalg.det(MA)

       ##5.定义单元应变矩阵并计算单元应变
       BE = 1 / 2 * A * np.array(
           [[b[0], 0, b[1], 0, b[2], 0], [0, c[0], 0, c[1], 0, c[2]], [c[0], b[0], c[1], b[1], c[2], b[2]]])
       be = np.dot(BE, ue)     #单元应变向量
       ##6.将应变添加到结构应变信息矩阵
       B[n - 1, :] = [n, be[0, 0], be[1, 0], be[2, 0]]
   return B
def solve_S(U, NODE, ELE, ELE_SUM, EL, PO):
   """
   #根据位移向量U求解结构应力
   返回值：结构应变信息矩阵S
   第一列：单元号 第二列~第四列：σx,σy,τxy
   """
   S = np.zeros((ELE_SUM, 4))  #初始化单元应力信息矩阵
   for n in range(1, ELE_SUM + 1):
       ##1.建立当前单元的位移向量
       ue = np.zeros((6, 1))  # 初始化单元位移向量ue
       for i in range(3):
           ue[2 * i, 0] = U[2 * ELE[n - 1, i] - 2]
           ue[2 * i + 1, 0] = U[2 * ELE[n - 1, i] - 1]
       ##2.整理当前单元节点坐标##
       x = np.array([0] * 6)
       y = np.array([0] * 6)  # 当前单元坐标数组
       for i in range(3):     # 将当前单元的坐标存入数组x,y，顺序为单元节点号矩阵内的顺序，编号从0开始(x0,x1,x2)
           x[i] = NODE[ELE[n - 1, i] - 1, 0]
           y[i] = NODE[ELE[n - 1, i] - 1, 1]
       x[3:6] = x[0:3]
       y[3:6] = y[0:3]        # 便于后续下标轮换调用
       ##3.计算插值函数系数
       a = np.array([0] * 3);
       b = np.array([0] * 3);
       c = np.array([0] * 3)  # a,b,c为插值函数系数，初始化
       for i in range(3):
           a[i] = x[i + 1] * y[i + 2] - x[i + 2] * y[i+1]  # ai=xjym-xmyj,下标i，j，m循环(a0,a1,a2)
           b[i] = y[i + 1] - y[i + 2]  # bi=yi-ym
           c[i] = -x[i + 1] + x[i + 2]  # ci=-xi+xm
       ##4.计算三角形单元面积A
       MA = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])
       A = 0.5 * np.linalg.det(MA)
       ##5.计算单元应力
       SE=np.zeros((3,6))  #初始化单元应力矩阵
       for i in range(3):
           SE[:,2*i:2*i+2]=EL/(2*(1-PO**2)*A)*np.array([ [b[i],PO*c[i]],[PO*b[i],c[i]],[(1-PO)/2*c[i],(1-PO)/2*b[i]]])
       se= np.dot(SE, ue)  #计算单元应力向量[σx,σy,τ,]T
       ##6.将应变添加到结构应变信息矩阵
       S[n - 1, :] = [n, se[0, 0], se[1, 0], se[2, 0]]
   return S



"""
   ###主函数部分###
"""
##一、输入参数
#————————————————————————————————————————
#1.材料参数
#1.材料参数
EL = 20  # 弹性模量
PO =0.1  # 泊松比
#2.几何信息
TK = 1  # 厚度
NODE_SUM = 7  # 节点数目
ELE_SUM = 5  # 单元数目
NODE = np.array([[0,1], [0,0], [2,0], [2, 1],[4,1],[4,0],[6,1]])  # 节点坐标，按节点号排列(节点号必须为整数)
ELE = np.array([[1,2,3], [3, 4, 1],[4,3,5],[5,6,4],[5,6,7]])  # 单元节点号，单元内逆时针排列
#3.荷载和边界条件
P_inp= np.array([[7,0,-1]])    #输入荷载信息，第一列为荷载作用的节点号，第二列为x方向荷载，第三列为y方向荷载
BC=np.array([[1,1,1],[2,1,1]])          #约束信息，第一列为节点号，二~三列分别为x和y向的约束情况，1为固定，0为自由
#4.平面问题类型
TYPE=2  #1为平面应变问题，其余值为平面应力问题

##二、调用函数
#————————————————————————————————————————
if TYPE==1:
   EL=EL/(1-PO**2)
   PO=PO/(1-PO)
#1.求节点位移
print('平面三角形单元有限元程序开始求解')
ASK=create_K( EL, PO, TK, NODE_SUM,ELE_SUM, NODE, ELE)  #生成结构刚度矩阵
print('生成刚度矩阵：成功')
P=create_P(P_inp,NODE_SUM)  #生成荷载矩阵
print('生成荷载矩阵：成功')
do_BC(BC)   #处理边界条件
U=s.sparse.linalg.spsolve(ASK,P)    #调用sparse.linalg.spsolve函数，解大型稀疏矩阵方程
print('求解节点位移：成功')
#2.由节点位移求应变、应力
Uout=solve_Uout(U,NODE_SUM)
B=solve_B(U,NODE,ELE,ELE_SUM)   #得到单元应变信息矩阵
S=solve_S(U, NODE, ELE, ELE_SUM, EL, PO)    #得到单元应力信息矩阵
##三、输出结果
#————————————————————————————————————————

fw =open('result.txt','w', encoding="utf-8")
fa=open('result.txt','a',encoding="utf-8")
fw.write("——————###结果文件###——————\n\n")
fa.write('##节点位移##\n节点号\tx\ty\n——————————————————\n')
np.savetxt(fa, Uout, '%.8f',delimiter='\t')
fa.write('\n\n##单元应变##\n单元号\tεx\tεy\tγxy\n——————————————————\n')
np.savetxt(fa, B, '%.8f',delimiter='\t')
fa.write('\n\n##单元应力##\n单元号\tσx\tσy\tτxy\n——————————————————\n')
np.savetxt(fa, S, '%.8f',delimiter='\t')