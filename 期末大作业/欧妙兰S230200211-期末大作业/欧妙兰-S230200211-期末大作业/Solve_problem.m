function [Displacement,Force]=Solve_problem(Stiffness_matrix,...
    Fixed_situation_x,Fixed_situation_y,External_load_x,External_load_y,Num_node)
%% 记录非固定节点以及方向
Record=[];%记录那一部分是固定的
for i=1:Num_node %判断是否固定x
    if Fixed_situation_x(i)==0 %只记录那些非固定的
        Record=[Record 2*i-1];
    end
    if Fixed_situation_y(i)==0
        Record=[Record 2*i];
    end
end
%% 外载组成列向量
External_load_sum=[];%外载矩阵组合成一起
for i = 1:2*Num_node
    if mod(i,2)==1
        External_load_sum(i)=External_load_x((i+1)/2);
    else
        External_load_sum(i)=External_load_y(i/2);
    end
end
%% 筛选出非固定节点的对应刚度矩阵
Free_External_load=[];%外载自由矩阵筛选
Free_stiffness=[];% 刚度矩阵筛选出自由的来
Num_record=size(Record,2);
for i =1:Num_record
    Free_External_load(i)=External_load_sum(Record(i));
    for j =1:Num_record
        Free_stiffness(i,j)=Stiffness_matrix(Record(i),Record(j));
    end
end
%% 求解出节点位移以及支反力
Free_Displacement=(Free_stiffness^-1)*Free_External_load';%求解得出自由位移
Displacement=zeros(2*Num_node,1);%整体位移 方便后续做和
for i =1:Num_record %重新组装为整体
    Displacement(Record(i))=Free_Displacement(i);
end
Force=Stiffness_matrix*Displacement;%得到各个节点所得到的支反力

end
    

