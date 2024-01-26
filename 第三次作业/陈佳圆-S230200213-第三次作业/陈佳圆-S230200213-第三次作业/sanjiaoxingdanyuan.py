import numpy as np
import scipy as s
'''
本次三角形单元的程序设计流程为：
①编写生成载荷向量的函数
②编写生成整体刚度矩阵的函数（为核心部分，先生成单元的刚度矩阵，再进行对号入座，组装整体的刚度矩阵）
③编写处理边界条件的函数（利用的是划0置1法）
④编写输出位移矩阵的函数
⑤编写输出B矩阵的函数
⑥编写输出应变矩阵的函数
⑦在主函数中指定各种材料参数与单元数目与节点位置信息等
'''
def create_P(P_inp,NODE_SUM):#由输入的节点荷载信息生成荷载向量P
   P=np.zeros((2*NODE_SUM,1))  #初始化荷载向量P
   k=1
   P2 = P_inp.astype(int)
   for i in P2[:,0]:
       P[2*i-2]=P_inp[k-1,1]#X方向的载荷
       P[2*i-1]=P_inp[k-1,2]#Y方向的载荷
       k=k+1
   return (P)
def create_K( TM, BSB, HD, NODE_SUM,ELE_SUM, NODE, ELE):  # 生成整体刚度矩阵
   ASK = np.zeros((2 * NODE_SUM, 2 * NODE_SUM))#生成一个初始的总体刚度矩阵，其行数和列数均为两倍节点数目，即自由度的数目
   for n in range(1,ELE_SUM+1):#整理当前单元节点坐标
       x = np.array([0] * 6);y = np.array([0] * 6)#生成当前的单元坐标数组，因为每个三角形单元均有三个节点，每个节点均有两个自由度
       for i in range(3):#将当前单元的坐标存入数组，顺序为单元节点号矩阵内的顺序，编号从0开始(x0,x1,x2)
           x[i] = NODE[ELE[n - 1, i] - 1, 0]
           y[i] = NODE[ELE[n - 1, i] - 1, 1]
       x[3:6] = x[0:3]
       y[3:6] = y[0:3]# 便于后续下标轮换调用
       MA=np.array([[1,x[0],y[0]],[1,x[1],y[1]],[1,x[2],y[2]]])#创建一个3×3的矩阵。第一列均为1
       A=0.5*np.linalg.det(MA)#计算三角形单元面积A
       a = np.array([0] * 3);#创建一个[0,0,0]的一维数组
       b = np.array([0] * 3);#创建一个[0,0,0]的一维数组
       c = np.array([0] * 3)#创建一个[0,0,0]的一维数组，即a,b,c为插值函数系数，对其进行初始化
       for i in range(3):#计算插值函数系数，利用拉格朗日插值
           a[i] = x[i + 1] * y[i + 2] - x[i + 2] * y[i+1]#ai=xjym-xmyj,下标i，j，m循环得到(a0,a1,a2)
           b[i] = y[i+1] - y[i + 2]#bi=yi-ym循环得到（b0,b1,b2)
           c[i] = -x[i+1] + x[i + 2]#ci=-xi+xm循环得到（c0,c1,c2)
       EK=np.zeros((6,6))#生成单元刚度矩阵
       for i in range(3):
           for j in range(3):
               k1=b[i]*b[j]+(1-BSB)/2*c[i]*c[j]
               k2=BSB*c[i]*b[j]+(1-BSB)/2*b[i]*c[j]
               k3=BSB*b[i]*c[j]+(1-BSB)/2*c[i]*b[j]
               k4=c[i]*c[j]+(1-BSB)/2*b[i]*b[j]
               EK[2*i:2*i+2,2*j:2*j+2]=TM*HD/(4*(1-BSB**2)*A)*np.array([[k1,k3],[k2,k4]])#生成各部分的值后组装成一个2×2的矩阵
       ask=np.zeros((2*NODE_SUM,2*NODE_SUM)) #临时总体刚度矩阵的初始化
       G=np.zeros((6,2*NODE_SUM))#定义转换矩阵G(6行，2倍节点数目列），并使其初始化
       for i in range(3):#将特定位置置为1，其余位置保持0
           G[2*i,2*ELE[n-1,i]-2]=1
           G[2*i+1,2*ELE[n -1,i]-1]=1
       GT=np.transpose(G)#对G矩阵进行转置
       ask=np.dot(np.dot(GT,EK),G)#对GT和EK矩阵进行乘法后再与G进行矩阵乘法
       ASK=ASK+ask#将单元刚度矩阵组装为整体刚度矩阵
   return ASK

def do_BT(BT):#处理边界条件，置1法对角元素
   global ASK#申明整体刚度矩阵为全局变量
   global P#申明载荷矩阵为全局变量
   k=1
   for i in BT[:,0]:#对载荷矩阵进行赋值
       if BT[k-1,1]==1:
           P[2*i-2]=0#X方向的载荷置为0
       if BT[k-1,2]==1:
           P[2*i-1]=0#Y方向的载荷置为0
       k=k+1
   k=1
   for i in BT[:, 0]:#对总体刚度矩阵进行赋值
       if BT[k-1, 1]==1:
           ASK[2*i-2,:]=0#将2*i-2整行置0
           ASK[:,2*i-2]=0#将2*i-2整列置0
           ASK[2*i-2,2*i-2]=1#将第2*i-2行和第2*i-2列的交叉位置元素设置为1
       if BT[k-1,2]==1:
           ASK[2*i-1,:] = 0#将2*i-1整行置0
           ASK[:,2*i-1] = 0#将2*i-1整列置0
           ASK[2*i-1,2*i-1]=1#将第2*i-1行和第2*i-1列的交叉位置元素设置为1
       k=k+1
def solve_WYout(WY,NODE_SUM):
   WYout=np.zeros((NODE_SUM,3))#生成初始位移矩阵
   for i in range(NODE_SUM):
       WYout[i,0]=i+1
       WYout[i,1]=WY[2*i]
       WYout[i, 2]=WY[2*i+1]
   return WYout#输出位移矩阵

def solve_B(WY,NODE,ELE,ELE_SUM):
   B = np.zeros((ELE_SUM, 4))#初始化单元应变信息矩阵
   for n in range(1, ELE_SUM + 1):#建立当前单元的位移向量
       wy = np.zeros((6, 1))#初始化单元位移向量wy
       for i in range(3):#生成单元位移向量
           wy[2 * i, 0] = WY[2 * ELE[n - 1, i] - 2]
           wy[2 * i + 1, 0] = WY[2 * ELE[n - 1, i] - 1]
       x = np.array([0] * 6)
       y = np.array([0] * 6)#当前单元坐标数组
       for i in range(3):#整理当前单元节点坐标
           x[i] = NODE[ELE[n - 1, i] - 1, 0]#将当前单元的坐标存入数组x，顺序为单元节点号矩阵内的顺序，编号从0开始(x0,x1,x2)
           y[i] = NODE[ELE[n - 1, i] - 1, 1]#将当前单元的坐标存入数组y，顺序为单元节点号矩阵内的顺序，编号从0开始(y0,y1,y2)
       x[3:6] = x[0:3]
       y[3:6] = y[0:3]#便于后续下标轮换调用
       a = np.array([0] * 3);
       b = np.array([0] * 3);
       c = np.array([0] * 3)#a,b,c为插值函数系数，初始化
       for i in range(3):#计算插值函数系数
           a[i] = x[i + 1] * y[i + 2] - x[i + 2] * y[i+1]#ai=xjym-xmyj,下标i，j，m循环(a0,a1,a2)
           b[i] = y[i + 1] - y[i + 2]#bi=yi-ym循环(b0,b1,b2)
           c[i] = -x[i + 1] + x[i + 2]#ci=-xi+xm(c0,c1,c2)
       MA = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])#组建计算三角形单元面积A的矩阵
       A = 0.5 * np.linalg.det(MA)#计算三角形单元面积A

       BE = 1 / 2 * A * np.array(
           [[b[0], 0, b[1], 0, b[2], 0], [0, c[0], 0, c[1], 0, c[2]], [c[0], b[0], c[1], b[1], c[2], b[2]]])#定义单元应变矩阵
       be = np.dot(BE, wy)#单元应变向量
       B[n - 1, :] = [n, be[0, 0], be[1, 0], be[2, 0]]#将应变添加到结构应变信息矩阵
   return B
def solve_S(WY, NODE, ELE, ELE_SUM, TM, BSB):
   S = np.zeros((ELE_SUM, 4))#初始化单元应力信息矩阵
   for n in range(1, ELE_SUM + 1):#建立当前单元的位移向量
       wy = np.zeros((6, 1))#初始化单元位移向量wy
       for i in range(3):#给单元的位移矩阵赋值
           wy[2 * i, 0] = WY[2 * ELE[n - 1, i] - 2]
           wy[2 * i + 1, 0] = WY[2 * ELE[n - 1, i] - 1]
       x = np.array([0] * 6)
       y = np.array([0] * 6)#整理当前单元节点坐标
       for i in range(3):#将当前单元的坐标存入数组x,y，顺序为单元节点号矩阵内的顺序，编号从0开始(x0,x1,x2)
           x[i] = NODE[ELE[n - 1, i] - 1, 0]
           y[i] = NODE[ELE[n - 1, i] - 1, 1]
       x[3:6] = x[0:3]
       y[3:6] = y[0:3]#便于后续下标轮换调用
       a = np.array([0] * 3);
       b = np.array([0] * 3);
       c = np.array([0] * 3)  # a,b,c为插值函数系数，初始化
       for i in range(3):#计算插值函数系数
           a[i] = x[i + 1] * y[i + 2] - x[i + 2] * y[i+1]#ai=xjym-xmyj,下标i，j，m循环(a0,a1,a2)
           b[i] = y[i + 1] - y[i + 2]#bi=yi-ym循环(b0,b1,b2)
           c[i] = -x[i + 1] + x[i + 2]#ci=-xi+xm循环(c0,c1,c2)
       MA = np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])#组建计算三角形单元面积A的矩阵
       A = 0.5 * np.linalg.det(MA)#计算三角形单元面积A
       SE=np.zeros((3,6))#初始化单元应力矩阵
       for i in range(3):
           SE[:,2*i:2*i+2]=TM/(2*(1-BSB**2)*A)*np.array([ [b[i],BSB*c[i]],[BSB*b[i],c[i]],[(1-BSB)/2*c[i],(1-BSB)/2*b[i]]])
       se= np.dot(SE, wy)#计算单元应力向量[σx,σy,τ,]T
       S[n - 1, :] = [n, se[0, 0], se[1, 0], se[2, 0]]#将应变添加到结构应变信息矩阵
   return S

TM = 1#弹性模量
BSB = 0.3#泊松比
HD = 1#厚度
NODE_SUM = 4#节点数目（从0开始）
ELE_SUM = 2#单元数目（从0开始）
NODE = np.array([[0,2], [0,0], [4,0], [4, 2],[8,2]])# 节点坐标，按节点号排列(节点号必须为整数)
ELE = np.array([[1,2,3], [3,4,1],[3,5,4]])# 单元节点号，单元内逆时针排列
P_inp= np.array([[1,0,-0.5],[2,0,-0.5]])#输入荷载信息，第一列为荷载作用的节点号，第二列为x方向荷载，第三列为y方向荷载
BT=np.array([[3,1,1],[4,1,1]])#约束信息，第一列为节点号，二~三列分别为x和y向的约束情况，1为固定，0为自由
TYPE=2#1为平面应变问题，其余值为平面应力问题
if TYPE==1:
   TM=TM/(1-BSB**2)
   BSB=BSB/(1-BSB)
P=create_P(P_inp,NODE_SUM)#生成荷载矩阵
ASK=create_K( TM, BSB, HD, NODE_SUM,ELE_SUM, NODE, ELE)#生成结构刚度矩阵
do_BT(BT)#处理边界条件
WY=s.sparse.linalg.spsolve(ASK,P)#调用sparse.linalg.spsolve函数，解大型稀疏矩阵方程
WYout=solve_WYout(WY,NODE_SUM)#得到节点位移矩阵
B=solve_B(WY,NODE,ELE,ELE_SUM)#得到单元应变信息矩阵
S=solve_S(WY, NODE, ELE, ELE_SUM, TM, BSB)#得到单元应力信息矩阵
print(B)
print(S)