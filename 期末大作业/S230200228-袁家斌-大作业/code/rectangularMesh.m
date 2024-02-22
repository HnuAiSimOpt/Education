
function [elemsID,elemNodes, nodeCoordinates] = ...
    rectangularMesh(lx, ly, num_x, num_y)

%% 此函数用于划分规则矩形单元
% elemsID  单元编码
% elemNodes 单个单元的四个节点编号为此数据的一行
% nodeCoodinates 每个节点的(编号,x,y)坐标为此数据的一行
% lx 矩形长
% ly 矩形宽
% num_x 矩形x方向单元数量
% num_y 矩形y方向单元数量

%% 初始化
nodeCoordinates = zeros((num_x+1)*(num_y+1), 3);
count_node = 1;

elemsID = num_x * num_y;
elemNodes = zeros(num_x*num_y, 4);
count_elem = 1;

%% 赋值
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

