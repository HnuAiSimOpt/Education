function [globalStiffness] = globalStiffness2D(GDof, elemNodes,...
    nodeCoordinates, thickness, constMatrix)

%% �˺�����������ȫ�ָնȾ���
% globalStiffness ȫ�ָնȾ���
% GDof ȫ�����ɶ�
% elemNodes ��Ԫ�ڵ���
% nodeCoordinates �ڵ�����
% thickness ��Ԫ���
% constMatrix ��������

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
    
    %���ӵ�Ԫ�նȾ���
    indexs = zeros(1, 6*length(index));
    indexs(1:6:end) = 6*index-5; indexs(2:6:end) = 6*index-4;
    indexs(3:6:end) = 6*index-3; indexs(4:6:end) = 6*index-2;
    indexs(5:6:end) = 6*index-1; indexs(6:6:end) = 6*index;
    
    globalStiffness(indexs, indexs) = ...
        globalStiffness(indexs, indexs) + elemStiff;
end
globalStiffness(6:6:end, 6:6:end) = eye(size(nodeCoordinates, 1));
end