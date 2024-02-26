function element_solid=Solve_element(overview)
for i=2986:2:6704
Element_Solid(i-2985,:)=strsplit(overview{i,1});  % 将元胞中每个元素的内容分来
end
Element_Solid(:,1)=[];
for i=2986:2:6704
    Element_solid(i-2985,1)=str2num(Element_Solid{i-2985,1});
    Element_solid(i-2985,2)=str2num(Element_Solid{i-2985,2});
end
for i=2987:2:6705
Element_Node(i-2986,:)=strsplit(overview{i,1});
end
Element_Node(:,1)=[];
for i=2987:2:6705
    Element_solid(i-2986,2)=str2num(Element_Node{i-2986,1});
    Element_solid(i-2986,3)=str2num(Element_Node{i-2986,2});
    Element_solid(i-2986,4)=str2num(Element_Node{i-2986,3});
    Element_solid(i-2986,5)=str2num(Element_Node{i-2986,4});
    Element_solid(i-2986,6)=str2num(Element_Node{i-2986,5});
    Element_solid(i-2986,7)=str2num(Element_Node{i-2986,6});
    Element_solid(i-2986,8)=str2num(Element_Node{i-2986,7});
    Element_solid(i-2986,9)=str2num(Element_Node{i-2986,8});
end
i=1;
for b=1:1860
    element_solid(b,1)=Element_solid(i,1);
    element_solid(b,2)=Element_solid(i,2);
    element_solid(b,3)=Element_solid(i,3);
    element_solid(b,4)=Element_solid(i,4);
    element_solid(b,5)=Element_solid(i,5);
    element_solid(b,6)=Element_solid(i,6);
    element_solid(b,7)=Element_solid(i,7);
    element_solid(b,8)=Element_solid(i,8);
    element_solid(b,9)=Element_solid(i,9);
    i=i+2;
end
for b=1:1800
    element_solid(b,1)=element_solid(b,1)-124;
end