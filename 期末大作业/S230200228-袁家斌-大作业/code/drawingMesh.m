function [] = drawingMesh(elemNodes, nodeCoordinates)

%% �˺������ڻ�������ͼ��

xy_coord = nodeCoordinates(:,2:end);
patch('Faces',elemNodes,'Vertices',xy_coord,'FaceColor', 'green');
axis equal; axis tight; axis off;
end

