function [] = drawingMesh(elemNodes, nodeCoordinates)

%% 此函数用于绘制网格图形

xy_coord = nodeCoordinates(:,2:end);
patch('Faces',elemNodes,'Vertices',xy_coord,'FaceColor', 'green');
axis equal; axis tight; axis off;
end

