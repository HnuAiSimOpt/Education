 function [] = postDisplacement(elemNodes, nodeCoordinates, disp, scale,...
    type, origin)

%% 此函数用于显示位移结果
% elemNodes 单元节点信息
% nodeCoordinates 节点坐标信息
% disp 节点位移信息
% type 位移分量信息 X 或 Y
% origin 是否显示变形前图， 是为1，否为0

xy_coord = nodeCoordinates(:,2:end);
%% 绘制变形后
switch type
    case 'X'
        d = disp(:,1);   type = 'X方向';
    case 'Y'
        d = disp(:,2);   type = 'Y方向';
    case 'Z'
        d = disp(:,3);   type = 'Z方向';
    case 'M'
        d = sqrt(disp(:,1).^2 + disp(:,2).^2 + disp(:,3).^2);
        type = '综合';
end

ticks = linspace(min(d), max(d), 11);
patch('Faces',elemNodes,'Vertices',xy_coord+disp*scale,...
     'FaceVertexCData', d, 'FaceColor', 'interp', 'FaceAlpha', 0.95);
colormap jet
colorbar('Ticks', ticks)

%% 绘制原模型
if origin == 1
    p = patch('Faces',elemNodes,'Vertices',xy_coord,...
         'FaceColor', [0 0 1], 'FaceAlpha', 0.1);
     p.LineStyle = ':';
end

%% 设置图形属性
axis equal; axis tight; axis off;hold on
title([type, '位移 （mm）'])
view(-37.5,30);      
end
