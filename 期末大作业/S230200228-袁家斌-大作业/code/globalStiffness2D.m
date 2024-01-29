function [globalStiffness] = globalStiffness2D(GDof, elemNodes,...
    nodeCoordinates, thickness, constMatrix)

%% 此函数用于生产全局刚度矩阵
% globalStiffness 全局刚度矩阵
% GDof 全局自由度
% elemNodes 单元节点编号
% nodeCoordinates 节点坐标
% thickness 单元厚度
% constMatrix 本构矩阵

%%
globalStiffness = sparse(GDof, GDof);
for e = 1:size(elemNodes, 1)
    nodes = elemNodes(e,:);
    index = zeros(1,length(nodes));
    for i = 1:length(nodes)
        index(i) = find(nodeCoordinates(:,1)==nodes(i));
    end
    coords = nodeCoordinates(index, 2:end);
    elemStiff = elemStiffness2D(coords, thickness, constMatrix);
    
    %叠加单元刚度矩阵
    indexs = zeros(1, 6*length(index));
    indexs(1:6:end) = 6*index-5; indexs(2:6:end) = 6*index-4;
    indexs(3:6:end) = 6*index-3; indexs(4:6:end) = 6*index-2;
    indexs(5:6:end) = 6*index-1; indexs(6:6:end) = 6*index;
    
    globalStiffness(indexs, indexs) = ...
        globalStiffness(indexs, indexs) + elemStiff;
end
globalStiffness(6:6:end, 6:6:end) = eye(size(nodeCoordinates, 1));
end