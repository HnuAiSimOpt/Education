
function [load] = distributeLoad(dofs, coords, elems, pressure )
%% 此函数用于生成分布式载荷
% load 载荷向量， dofs * 1
% dofs 自由度数， 节点数 * 每个节点的自由度数
% coords 节点坐标矩阵 节点数 * 3  第1列为节点编号，2、3列为x,y坐标
% elems 单元节点集合

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