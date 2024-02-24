clear
clc
%% 读取题目信息数据
Input
%% 单元刚度矩阵的构建与全局刚度矩阵的装配
Dimension_Matrix=2*Num_Node;%矩阵维数
K_global=zeros(Dimension_Matrix);%初始化全局矩阵
Information_element=[];%接受单元信息 长度 角度值 且按单元顺序
for i = 1:Num_Element%全局矩阵装配
    [stiffness_element,length_element,c,s]=Element_stiffness_matrix(ID_Element(i,4),ID_Element(i,5),...
        ID_Node(ID_Element(i,2),2),ID_Node(ID_Element(i,2),3),...
        ID_Node(ID_Element(i,3),2),ID_Node(ID_Element(i,3),3));
    K_global=K_global+Assembly_stiffness_matrix(stiffness_element,...
        ID_Element(i,2),ID_Element(i,3),Num_Node);
    Information_element=[Information_element;length_element,c,s];
end
%% 求解桁架整体位移以及单元固定支反力
[Displacement,Force]=Solve_problem(K_global,ID_Node(:,6),ID_Node(:,7),ID_Node(:,4),ID_Node(:,5),Num_Node);
%% 求解题目要求竖直位移
Displacement_y=[];
Displacement_x=[];
for i =1:Dimension_Matrix %竖直位移
    if mod(i,2)==0
        Displacement_y=[Displacement_y;Displacement(i)];
    else
        Displacement_x=[Displacement_x;Displacement(i)];
    end
end
[max_Displacement_y,max_Displacement_pos]=min(Displacement_y);
%% 后处理求解单元应力和应变
[Stress,Strain]=Post_treatment(Displacement,ID_Element,Information_element);
%% 求解题目最大拉压应力
[max_stress,max_pos]=max(Stress);
[min_stress,min_pos]=min(Stress);
%% 输出图像结果
Output
