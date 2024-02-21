 function [] = postDisplacement(elemNodes, nodeCoordinates, disp, scale,...
    type, origin)

%% �˺���������ʾλ�ƽ��
% elemNodes ��Ԫ�ڵ���Ϣ
% nodeCoordinates �ڵ�������Ϣ
% disp �ڵ�λ����Ϣ
% type λ�Ʒ�����Ϣ X �� Y
% origin �Ƿ���ʾ����ǰͼ�� ��Ϊ1����Ϊ0

xy_coord = nodeCoordinates(:,2:end);
%% ���Ʊ��κ�
switch type
    case 'X'
        d = disp(:,1);   type = 'X����';
    case 'Y'
        d = disp(:,2);   type = 'Y����';
    case 'Z'
        d = disp(:,3);   type = 'Z����';
    case 'M'
        d = sqrt(disp(:,1).^2 + disp(:,2).^2 + disp(:,3).^2);
        type = '�ۺ�';
end

ticks = linspace(min(d), max(d), 11);
patch('Faces',elemNodes,'Vertices',xy_coord+disp*scale,...
     'FaceVertexCData', d, 'FaceColor', 'interp', 'FaceAlpha', 0.95);
colormap jet
colorbar('Ticks', ticks)

%% ����ԭģ��
if origin == 1
    p = patch('Faces',elemNodes,'Vertices',xy_coord,...
         'FaceColor', [0 0 1], 'FaceAlpha', 0.1);
     p.LineStyle = ':';
end

%% ����ͼ������
axis equal; axis tight; axis off;hold on
title([type, 'λ�� ��mm��'])
view(-37.5,30);      
end
