%%
%解abaqus三角形单元问题
%主程序
clear
clc

elements = load("main_3_elements.txt");
cast_nodes = 1:334:11023;
load_nodes = 334:334:11356;
nodes = load("main_3_nodes.txt");

[nodes_num, ~] = size(nodes);
[element_num, ~] = size(elements);
[~, load_nodes_num] = size(load_nodes);
[~, cast_nodes_num] = size(cast_nodes);
%%
E  = 210000;
NU = 0.3;
t  =1;
total_force_x = 100;
total_force_y = 10;
average_force_x = total_force_x / load_nodes_num;
average_force_y = total_force_y / load_nodes_num;
%%
format long
%%定义刚度矩阵
K   = zeros(nodes_num*2, nodes_num*2);
Patch_xy = zeros(6, element_num);

%%获取刚度矩阵K
for i = 1:element_num
    %%获取单元刚度矩阵
    node_1_index = elements(i, 2);
    node_2_index = elements(i, 3);
    node_3_index = elements(i, 4);
    node_1_coordinate = nodes(node_1_index, 2:3);
    node_2_coordinate = nodes(node_2_index, 2:3);
    node_3_coordinate = nodes(node_3_index, 2:3); 
    Patch_xy(:,i) = [node_1_coordinate(1); node_2_coordinate(1); node_3_coordinate(1);node_1_coordinate(2); node_2_coordinate(2); node_3_coordinate(2)];
    k_element = LinearTriangleElementStiffness(E,NU,t,...
                                              node_1_coordinate(1),node_1_coordinate(2),...
                                              node_2_coordinate(1),node_2_coordinate(2),...
                                              node_3_coordinate(1),node_3_coordinate(2),...
                                              1);
    %%组装刚度矩阵
    disp("组装"+string(i)+"刚度矩阵");
    nodes_element = [node_1_index, node_2_index, node_3_index];
    %可以调用LinearTriangleAssemble但是效率很低
    for row = 1:3
        for col = 1:3
            K((nodes_element(row)-1)*2+1:nodes_element(row)*2, (nodes_element(col)-1)*2+1:nodes_element(col)*2) = ...
                K((nodes_element(row)-1)*2+1:nodes_element(row)*2, (nodes_element(col)-1)*2+1:nodes_element(col)*2)+k_element((row-1)*2+1:row*2,(col-1)*2+1:col*2);
        end
    end
end
%%获取载荷矩阵F
%%
F   = zeros(nodes_num*2, 1);
KF  = zeros(nodes_num*2, nodes_num*2+1);
F(load_nodes*2, 1) = average_force_y;
F(load_nodes*2-1, 1) = average_force_x;
%%获得增广矩阵
KF(:,1:nodes_num*2) = K;
KF(:,nodes_num*2+1) = F;
%%施加边界条件，获得施加边界U后的KF
del_index = zeros(1, cast_nodes_num*2);
del_index(1:cast_nodes_num) = cast_nodes*2;
del_index(cast_nodes_num+1:cast_nodes_num*2) = cast_nodes*2 - 1;
KF(del_index,:) = [];
KF(:,del_index) = [];

%%
[KUF_num, ~] = size(KF);
U = KF(:, 1:KUF_num)\KF(:,KUF_num+1);

for i = 1:cast_nodes_num
    cast_index = cast_nodes(i);
    U = [U(1:(cast_index-1)*2,:); 0; 0; U(((cast_index-1)*2+1):end,:)];
end

% scatter(nodes(:,2), nodes(:,3),50,U(1:2:end),"filled");
% colorbar
%%
s = zeros(3,element_num);

for i =  1:element_num
    node_1_index = elements(i, 2);
    node_2_index = elements(i, 3);
    node_3_index = elements(i, 4);
    u = [U(node_1_index*2-1);U(node_1_index*2);U(node_2_index*2-1);U(node_2_index*2);U(node_3_index*2-1);U(node_3_index*2);];
    s(:,i) = LinearTriangleElementPStresses(LinearTriangleElementStresses(E, NU, Patch_xy(1,i), Patch_xy(4,i), Patch_xy(2,i), Patch_xy(5,i),Patch_xy(3,i), Patch_xy(6,i), 1, u),2);
end
%%
clf
hold on

patch(Patch_xy(1:3,:), Patch_xy(4:6,:),s(1,:));
colorbar

