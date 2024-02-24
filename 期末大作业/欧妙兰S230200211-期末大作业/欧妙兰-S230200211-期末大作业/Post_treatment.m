function [Stress,Strain]=Post_treatment(Displacement,ID_Element,Information_element)
%% 利用公式将全局位移转化为局部位移
Num_element=size(ID_Element,1);
Strain=[];%初始化单元局部位移
for i =1:Num_element
    u=[Information_element(i,2) Information_element(i,3) 0 0;...
        0 0 Information_element(i,2) Information_element(i,3)]...
        *[Displacement(2*ID_Element(i,2)-1);Displacement(2*ID_Element(i,2));...
        Displacement(2*ID_Element(i,3)-1);Displacement(2*ID_Element(i,3))];
    Strain=[Strain;(u(2)-u(1))/Information_element(i,1)];
end
Stress=ID_Element(:,4).*Strain;%直接运用公式求应力

end

