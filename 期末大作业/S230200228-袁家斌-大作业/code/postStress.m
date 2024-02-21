
function [] = postStress(elemNodes, nodeCoordinates, scale, disp, ...
    stress, type, h)
%% �˺���������ʾλ�ƽ��
% elemNodes ��Ԫ�ڵ���Ϣ
% nodeCoordinates �ڵ�������Ϣ
% disp λ�Ʒ���
% stress �ڵ�Ӧ����Ϣ
% type Ӧ��������sigma_x, sigma_y, tao_xy, mises
% �غ�ȷ���Ӧ��λ��

xy_coord = nodeCoordinates(:,2:end);
%% ���Ʊ��κ�
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

%% ����ͼ������
axis equal; axis tight; axis off;hold on
title(['�ں��',num2str(h),'mm����', 'Ӧ�� ��MPa��'])
view(-37.5,30);     
end

