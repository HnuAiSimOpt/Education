function z = Quad2D4Node_Assembly(KK,k,i,j,m,n)
%该函数进行单元刚度矩阵的组装
%输入单元刚度矩阵 k，单元的节点编号 i、j、m、n
%输出整体刚度矩阵 KK

DOF(1)=2*i-1;
DOF(2)=2*i;
DOF(3)=2*j-1;
DOF(4)=2*j;
DOF(5)=2*m-1;
DOF(6)=2*m;
DOF(7)=2*n-1;
DOF(8)=2*n;
for n1=1:8
    for n2=1:8
        KK(DOF(n1),DOF(n2))= KK(DOF(n1),DOF(n2))+k(n1,n2);
    end
end
z=KK;

