function [K] = funct_Stiff_Beam3D(nP,nE,eNodes,xyz,...
    A,E,ks,G,Iy,Iz,kt,direcVec)
%%%% 计算梁结构刚度矩阵
nDof = nP*6;
K = sparse(nDof,nDof);

gaussLoc = 1/sqrt(3)*[-1; 1];
gaussWts = [1; 1];

for i = 1:nE
    eNodei = eNodes(i,:); %%%%% eNodei = [node 1st, node 2nd]
    nnode = length(eNodei);%%%% 每个元素的节点数 = 2
    
    %%%% det(J) and J^-1
    [detJ, invJ] = detJ_invJ(xyz,eNodei);
    
    %%%% 简单计算剪切和轴向刚度矩阵
    %%%% 组装
    [keLoc] = stiff_shear_axial(invJ,nnode,A,E,detJ,ks,G,i);
    
    %%%% 完整计算弯曲和扭转刚度矩阵
    %%%% 组装
    [keLoc2] = stiff_bend_tors(keLoc,gaussLoc,gaussWts,...
        detJ,invJ,nnode,kt,G,E,Iy,Iz,i);
    %%%%% 局部刚度矩阵计算完成
    %%%%% 元素，转移到全局坐标系
    
    %%%%%% 刚度的变换矩阵
    [keG,eDofs] = transform_stiff(direcVec,keLoc2,i,eNodei,nP);
    
    %%%%% 全局刚度矩阵
    K(eDofs,eDofs) = K(eDofs,eDofs) + keG;
end
end

