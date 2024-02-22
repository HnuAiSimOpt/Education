%%
%载入节点
clear
clc
format short
Nodes = load("Nodes.txt");
Load_nodes = load("Load_Nodes.txt")';
Boundary_nodes = load("Boundary_Nodes.txt")';
C3D4_Elements = load("3D4C_Elements.txt");

[Total_C3D4, ~] = size(C3D4_Elements);
Nodes_num = length(Nodes);
Dof = 3;
Total_Force = [0 ,-100, 0];


Total_Dof = Dof*Nodes_num;
%%
figure1 = figure();
hold on;
scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4),"g.");
axis equal;

%边界节点
Boundary_nodes_coor = Nodes(Boundary_nodes,2:4);
scatter3(Boundary_nodes_coor(:,1),Boundary_nodes_coor(:,2),Boundary_nodes_coor(:,3),"r*");
%加载节点
Load_nodes_coor = Nodes(Load_nodes,2:4);
scatter3(Load_nodes_coor(:,1),Load_nodes_coor(:,2),Load_nodes_coor(:,3),"yo");
view([45,45])

%%
%材料属性编辑
elastic = 210000;
poisson = 0.3;
iopt = 4;
D = Element_D(iopt,elastic,poisson);
%%
K = zeros(Total_Dof,Total_Dof);

%%
%四面体单元
for i = 1:Total_C3D4
    element_nodes_index = C3D4_Elements(i,2:5);
    four_nodes_matrix = Nodes(element_nodes_index,2:end);
    disp("组装四面体单元"+string(i)+"个刚度矩阵");
    k = C3D4_K(D, four_nodes_matrix);
    for row = 1:4
        row_index = element_nodes_index(row);%刚度矩阵行节点编号
        for col = 1:4
            col_index = element_nodes_index(col);%刚度矩阵列节点编号
            K((3*row_index-2):row_index*3,(3*col_index-2):col_index*3) = ...
                K((3*row_index-2):row_index*3,(3*col_index-2):col_index*3)...
                + k((3*row-2):3*row,(3*col-2):3*col);
        end
    end
end
%%
F = Load_Apply(Load_nodes, Nodes_num, Dof, Total_Force);
Constrain_dofs = zeros(length(Boundary_nodes),3);
for i = 1:Dof
    Constrain_dofs(:,i) = (Boundary_nodes-1)*Dof+i;
end
switch Dof
    case 1
        Constrain = Constrain_dofs(:,1);
    case 2
        Constrain = [Constrain_dofs(:,1);Constrain_dofs(:,2)];
    case 3
        Constrain = [Constrain_dofs(:,1);Constrain_dofs(:,2);Constrain_dofs(:,3)];
end
K_constrain = K;
F_constrain = F;
K_constrain(Constrain,:) = [];
F_constrain(Constrain,:) = [];
K_constrain(:,Constrain) = [];

%%
%求位移节点矩阵
disp("开始计算");
U_ = K_constrain\F_constrain;
%%
%返回节点位移
U = U_;
for i = 1:length(Boundary_nodes)
    index = Boundary_nodes(i);
    forward_ = U(1:(index-1)*3,:);
    backward_ = U((index-1)*3+1:end,:); 
    U = [forward_;0;0;0;backward_];
end
%%
U1 = U(1:3:end);
U2 = U(2:3:end);
U3 = U(3:3:end);
figure2 = figure();
scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4),20,U1,'fill');
figure3 = figure();
scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4),20,U2,'fill');
figure4 = figure();
scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4),20,U3,'fill');
%%
%保存计算文件
fileID = fopen("U1.txt", 'w');
fprintf(fileID, '%f\n',U1);
fclose(fileID);
fileID = fopen("U2.txt", 'w');
fprintf(fileID, '%f\n',U2);
fclose(fileID);
fileID = fopen("U3.txt", 'w');
fprintf(fileID, '%f\n',U3);
fclose(fileID);
%%
tetramesh(C3D4_Elements(:,2:5), Nodes(:,2:4),'FaceColor','white');
hold on
%%
figure3 = figure();
scatter3(Nodes(:,2),Nodes(:,3),Nodes(:,4),20,U3,'fill');








