
function [load] = distributeLoad(dofs, coords, elems, pressure )
%% �˺����������ɷֲ�ʽ�غ�
% load �غ������� dofs * 1
% dofs ���ɶ����� �ڵ��� * ÿ���ڵ�����ɶ���
% coords �ڵ�������� �ڵ��� * 3  ��1��Ϊ�ڵ��ţ�2��3��Ϊx,y����
% elems ��Ԫ�ڵ㼯��

%%
[weights, points] = gaussIntegration(2);
load = sparse(dofs, 1);

for e = 1:size(elems, 1)
    nodes = elems(e,:);
    index = zeros(1,length(nodes));
    for i = 1:length(nodes)
        index(i) = find(coords(:, 1)==nodes(i));
    end
    elemcoords = coords(index, 2:3);
    
    for g = 1:size(weights)
    	xi = points(g, 1);
        eta = points(g, 2);
        [N, dN] = shapeFun2D(xi, eta, 'Q4');
        [J, ~] = jacobian2D(elemcoords, dN);
        indexs = 6 * index - 3; 
        load(indexs) = load(indexs) + N * pressure * det(J) * weights(i);
    end
end
end