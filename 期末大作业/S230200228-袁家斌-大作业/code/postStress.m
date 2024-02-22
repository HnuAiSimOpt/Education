
function [] = postStress(elemNodes, nodeCoordinates, scale, disp, ...
    stress, type, h)
%% 此函数用于显示位移结果
% elemNodes 单元节点信息
% nodeCoordinates 节点坐标信息
% disp 位移分量
% stress 节点应力信息
% type 应力分量，sigma_x, sigma_y, tao_xy, mises
% 沿厚度方向应力位置

xy_coord = nodeCoordinates(:,2:end);
%% 绘制变形后
switch type
    case 'S11'
        s = stress(:,1);   type = 'S11';
    case 'S22'
        s = stress(:,2);   type = 'S22';
    case 'S12'
        s = stress(:,6);   type = 'S12';
    case 'S13'
        s = stress(:,4);   type = 'S13';
    case 'S23'
        s = stress(:,5);   type = 'S23';
    case 'M'
        s = sqrt(stress(:,1).^2 + stress(:,2).^2 - stress(:,1).*stress(:,2)...
        + 3*stress(:,6).^2);
        type = 'Mises';
end

ticks = linspace(min(s), max(s), 11);
patch('Faces',elemNodes,'Vertices',xy_coord+disp*scale,...
     'FaceVertexCData', s, 'FaceColor', 'interp', 'FaceAlpha', 0.95);
colormap jet
colorbar('Ticks', ticks)

%% 设置图形属性
axis equal; axis tight; axis off;hold on
title(['在厚度',num2str(h),'mm处的', '应力 （MPa）'])
view(-37.5,30);     
end

