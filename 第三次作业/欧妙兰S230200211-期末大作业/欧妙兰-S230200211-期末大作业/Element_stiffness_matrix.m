function [k_global,Length_element,c,s] = Element_stiffness_matrix(E,A,x1,y1,x2,y2)
%% 根据公式求解一个单元的全局刚度矩阵
Length_element=((y2-y1)^2+(x2-x1)^2)^0.5;
c=(x2-x1)/Length_element;
s=(y2-y1)/Length_element;
Transformation_matrix=[c^2 c*s -c^2 -s*c;s*c s^2 -s*c -s^2;-c^2 -s*c c^2 s*c;-s*c -s^2 s*c s^2];
k_element=A*E/Length_element;
k_global=k_element*Transformation_matrix;
end

