%%%%一阶六面体单元线弹性材料应力-应变矩阵D
function [D]=LinearIsotropicD(E,u) 
D=E/((1+u)*(1-2*u))*[1-u u u 0 0 0;u 1-u u 0 0 0;u u 1-u 0 0 0;0 0 0 (1-2*u)/2 0 0;0 0 0 0 (1-2*u)/2 0;0 0 0 0 0 (1-2*u)/2];
end