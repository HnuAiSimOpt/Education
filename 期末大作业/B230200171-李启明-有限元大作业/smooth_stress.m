function [Node1,Element1,stress1]= smooth_stress(stress,node,element)
%stress:单元高斯应力点应力值
%node:节点坐标
%element：单元连接矩阵

Node1 = [];
Element1 = [];
stress1 = [];
count = 0;

%单元重排列
for i = [1:size(element,1)]
    index_i = element(i,:);
    if sum(isnan(index_i))==0
        Node1 = [Node1; node(index_i,:)];
        Element1 = [Element1; count+1:count+8];
        stress1 = [stress1; ones(8,1)*stress(i)];
        count = count+8;
    end
end

stress_max = max(stress);
stress_min = min(stress);

Node_real=node; %独立坐标
% %容差
tol = 1e-6;
for i = [1:size(Node_real,1)]
    for j = [i+1:size(Node_real,1)]
        if abs(Node_real(j,1)-Node_real(i,1))<tol&&abs(Node_real(j,2)-Node_real(i,2))<tol&&abs(Node_real(j,3)-Node_real(i,3))<tol
             Node_real(j,:) = [inf,inf,inf]; 
        end
    end
end
Node_real(Node_real(:,1)==inf,:) = [];
Node_real = sortrows(Node_real);

%对每一个独立坐标寻找有多少个应力值
for i = [1:size(Node_real,1)]
    index = find(abs(Node_real(i,1)-Node1(:,1))<tol& abs(Node_real(i,2)-Node1(:,2))<tol&abs(Node_real(i,3)-Node1(:,3))<tol);
    stress_i = stress1(index);
    delta = (max(stress_i)-min(stress_i))/(stress_max-stress_min);
    if delta < 1
        stress1(index) = mean(stress_i);
    end
    
end
end