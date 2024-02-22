function [element_stiffness_matrix,B,D] = TetrahedronElementStiffiness(node_coordinate,element_node,E,mu)
[element_number,n_ele]=size(element_node);

dN_dlsf=[1	0	0	-1;
       0	1	0	-1;
       0	0	1	-1];
%求D矩阵
D=E*(1-mu)/(1+mu)/(1-2*mu)*[1 mu/(1-mu) mu/(1-mu) 0 0 0;
                            mu/(1-mu) 1 mu/(1-mu) 0 0 0;
                            mu/(1-mu) mu/(1-mu) 1 0 0 0;
                            0 0 0 (1-2*mu)/2/(1-mu) 0 0;
                            0 0 0 0 (1-2*mu)/2/(1-mu) 0;
                            0 0 0 0 0 (1-2*mu)/2/(1-mu)];
element_stiffness_matrix=zeros(3*n_ele,3*n_ele,element_number);
B=zeros(6,3*n_ele,element_number);
for i=1:element_number
    %% 求解B矩阵
        %组装为行列式
        node=node_coordinate(element_node(i,:),:);
        %求相关常数
        dN_dl=dN_dlsf;
        J=dN_dl*node;
        dN_dx=J\dN_dl;
        %组装B矩阵  
        for j=1:n_ele
            B(:,3*j-2:3*j,i)=[dN_dx(1,j) 0 0;
                      0 dN_dx(2,j) 0;
                      0 0 dN_dx(3,j);
                      dN_dx(2,j) dN_dx(1,j) 0;
                      0 dN_dx(3,j) dN_dx(2,j);
                      dN_dx(3,j) 0 dN_dx(1,j)];
        end
        %求解单元刚度矩阵
        element_stiffness_matrix(:,:,i)=B(:,:,i)'*D*B(:,:,i)*abs(det(J))/6;
end