% 绘制初始节点位置
figure;
plot(nodes(:, 1), nodes(:, 2), 'bo', 'MarkerSize', 4);
hold on;

% 绘制位移后的节点位置
plot(nodes_new(:, 1), nodes_new(:, 2), 'rx', 'MarkerSize', 4);

% 绘制单元连接关系
for i = 1:size(elements, 1)
    element_nodes = elements(i, :);
    %plot([nodes(element_nodes(1), 1), nodes(element_nodes(2), 1)], ...
        %[nodes(element_nodes(1), 2), nodes(element_nodes(2), 2)], 'b-');
        plot([nodes(element_nodes(1), 1), nodes(element_nodes(2), 1), nodes(element_nodes(3), 1)], ...
     [nodes(element_nodes(1), 2), nodes(element_nodes(2), 2), nodes(element_nodes(3), 2)], 'b-');
    plot([nodes_new(element_nodes(1), 1), nodes_new(element_nodes(2), 1), nodes_new(element_nodes(3), 1)], ...
     [nodes_new(element_nodes(1), 2), nodes_new(element_nodes(2), 2), nodes_new(element_nodes(3), 2)], 'r-');
end

xlabel('X 坐标');
ylabel('Y 坐标');
title('位移前后单元位置对比');
legend('初始位置', '位移后位置');
grid on;
