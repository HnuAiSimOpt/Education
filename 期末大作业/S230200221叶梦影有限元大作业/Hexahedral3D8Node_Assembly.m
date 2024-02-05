function z = Hexahedral3D8Node_Assembly(KK,k,i,j,l,m,n,o,p,q)
%该函数进行单元刚度矩阵的组装
%输入单元刚度矩阵 k
%输入单元的节点编号 i、j、l、m、n、o、p、q
%输出整体刚度矩阵 KK

DOF=[3*i-2:3*i,3*j-2:3*j,3*l-2:3*l,3*m-2:3*m,3*n-2:3*n,3*o-2:3*o,3*p-2:3*p,3*q-2:3*q];
for n1=1:24
    for n2=1:24
        KK(DOF(n1),DOF(n2))=KK(DOF(n1),DOF(n2))+k(n1,n2);
    end
end
z=KK;