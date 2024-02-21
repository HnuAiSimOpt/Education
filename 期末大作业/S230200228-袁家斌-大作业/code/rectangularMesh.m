
function [elemsID,elemNodes, nodeCoordinates] = ...
    rectangularMesh(lx, ly, num_x, num_y)

%% �˺������ڻ��ֹ�����ε�Ԫ
% elemsID  ��Ԫ����
% elemNodes ������Ԫ���ĸ��ڵ���Ϊ�����ݵ�һ��
% nodeCoodinates ÿ���ڵ��(���,x,y)����Ϊ�����ݵ�һ��
% lx ���γ�
% ly ���ο�
% num_x ����x����Ԫ����
% num_y ����y����Ԫ����

%% ��ʼ��
nodeCoordinates = zeros((num_x+1)*(num_y+1), 3);
count_node = 1;

elemsID = num_x * num_y;
elemNodes = zeros(num_x*num_y, 4);
count_elem = 1;

%% ��ֵ
for y = linspace(0, ly, num_y+1)
    for x = linspace(0, lx, num_x+1)
        nodeCoordinates(count_node,:) = [count_node, x, y];
        if ~(x==lx || y==ly)
            elemNodes(count_elem,:) = [count_node, count_node+1,...
                count_node+2+num_x, count_node+1+num_x];
            count_elem = count_elem + 1;
        else
        end
        count_node = count_node + 1;
    end
end

