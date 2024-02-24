function [Global_stiffness_matrix]=Assembly_stiffness_matrix(Element_stiffness_matrix,I,J,Num_node)
%% 刚度矩阵的装配 需要把单个全局刚度矩阵扩充到整体全局刚度矩阵之中 须在循环中依次进行
Global_stiffness_matrix=zeros(2*Num_node);%初始化刚度矩阵
Identifier=[2*I-1,2*I,2*J-1,2*J];
for i=1:4%通过两次循环加入到全局刚度矩阵之中
    for j=1:4
        Global_stiffness_matrix(Identifier(i),Identifier(j))=Element_stiffness_matrix(i,j);
    end
end
end
