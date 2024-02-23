function [direcVec] = directVectors(nE,eNodes,xyz,refP)
%%%% 计算所有元素的方向向量
direcVec = zeros(nE,9);
for i = 1:nE
    node1 = eNodes(i,1);%%%% 开始节点
    node2 = eNodes(i,2);%%%% 终止节点
    dxyz = xyz(node2,:) - xyz(node1,:);
    [Le] = lenVect(dxyz);
    
    Vri = 1/Le*dxyz;%%%% 单位向量
    Vs0 = refP - xyz(node1,:);
    [le_Vs0] = lenVect(Vs0);
    Vs0 = 1/le_Vs0*Vs0;%%%% 单位向量
    
    Vti = cross(Vri,Vs0);
    [le_Vti] = lenVect(Vti);
    Vti = 1/le_Vti*Vti;%%%% 单位向量
    
    Vsi = cross(Vti,Vri);
    [le_Vsi] = lenVect(Vsi);
    Vsi = 1/le_Vsi*Vsi;%%%% 单位向量
    
    Vri = cross(Vsi,Vti);
    [le_Vri] = lenVect(Vri);
    Vri = 1/le_Vri*Vri;%%%% 单位向量
    
    
    direcVec(i,:) = [Vri,Vsi,Vti];
end
end

