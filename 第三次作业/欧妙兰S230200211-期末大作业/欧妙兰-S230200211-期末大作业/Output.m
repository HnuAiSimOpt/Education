%%绘制单元及节点编号
% subplot(2,2,1);
title('单元及节点编号图');
for i=1:Num_Element
    axis([-1 110 0 80]);
    line([ID_Node(ID_Element(i,2),2) ID_Node(ID_Element(i,3),2)],...
        [ID_Node(ID_Element(i,2),3) ID_Node(ID_Element(i,3),3)],'color','k');
    %%绘制变形前单元
    text((ID_Node(ID_Element(i,2),2)+ID_Node(ID_Element(i,3),2))/2,...
        (ID_Node(ID_Element(i,2),3)+ID_Node(ID_Element(i,3),3))/2,...
        ['E',num2str(i)],'color','r','FontSize',7)
    hold on
end

for j=1:Num_Node
    %%绘制变形前节点
    plot(ID_Node(j,2),ID_Node(j,3),'k.','MarkerSize',15);
    text(ID_Node(j,2),ID_Node(j,3),['N',num2str(j)],'color','b','FontSize',7)
end

%%绘制变形图
% subplot(2,2,2);
%title('结构原图与变形图');
%for i=1:Num_Element
    %axis([-1 110 0 80]);
     %%绘制变形前单元
    %line([ID_Node(ID_Element(i,2),2) ID_Node(ID_Element(i,3),2)],...
    %    [ID_Node(ID_Element(i,2),3) ID_Node(ID_Element(i,3),3)],'color','k');
    %hold on
    %%绘制变形后单元
   % line([ID_Node(ID_Element(i,2),2)+Displacement(2*ID_Element(i,2)-1) ...
        %ID_Node(ID_Element(i,3),2)+Displacement(2*ID_Element(i,3)-1)],...
       % [ID_Node(ID_Element(i,2),3)+Displacement(2*ID_Element(i,2)) ...
        %ID_Node(ID_Element(i,3),3)+Displacement(2*ID_Element(i,3))],'color','r');
%end

%for j=1:Num_Node
    %%绘制变形前节点
    %plot(ID_Node(j,2),ID_Node(j,3),'k.','MarkerSize',15);
     %%绘制变形后节点
    %plot(ID_Node(j,2)+Displacement(2*j-1),ID_Node(j,3)+Displacement(2*j),'k.','MarkerSize',15);
%end

%%绘制最大拉压应力图
%subplot(2,2,3);
%title('最大拉/压应力图');
%for i=1:Num_Element
    %axis([-1 110 0 80]);
    %%绘制变形后单元
    %line([ID_Node(ID_Element(i,2),2)+Displacement(2*ID_Element(i,2)-1) ...
       % ID_Node(ID_Element(i,3),2)+Displacement(2*ID_Element(i,3)-1)],...
        %[ID_Node(ID_Element(i,2),3)+Displacement(2*ID_Element(i,2)) ...
       % ID_Node(ID_Element(i,3),3)+Displacement(2*ID_Element(i,3))],'color','r');
    %hold on
%end

%for j=1:Num_Node
     %%绘制变形后节点
   % plot(ID_Node(j,2)+Displacement(2*j-1),ID_Node(j,3)+Displacement(2*j),'k.','MarkerSize',15);
%end
%绘制最大拉应力
%line([ID_Node(ID_Element(max_pos,2),2)+Displacement(2*ID_Element(max_pos,2)-1) ...
       % ID_Node(ID_Element(max_pos,3),2)+Displacement(2*ID_Element(max_pos,3)-1)],...
       % [ID_Node(ID_Element(max_pos,2),3)+Displacement(2*ID_Element(max_pos,2)) ...
       % ID_Node(ID_Element(max_pos,3),3)+Displacement(2*ID_Element(max_pos,3))],'color','r','linewidth',5);
%绘制最大压应力
%line([ID_Node(ID_Element(min_pos,2),2)+Displacement(2*ID_Element(min_pos,2)-1) ...
       % ID_Node(ID_Element(min_pos,3),2)+Displacement(2*ID_Element(min_pos,3)-1)],...
        %[ID_Node(ID_Element(min_pos,2),3)+Displacement(2*ID_Element(min_pos,2)) ...
        %ID_Node(ID_Element(min_pos,3),3)+Displacement(2*ID_Element(min_pos,3))],'color','b','linewidth',5);

%text((ID_Node(ID_Element(max_pos,2),2)+Displacement(2*ID_Element(max_pos,2)-1)+ ...
      %  ID_Node(ID_Element(max_pos,3),2)+Displacement(2*ID_Element(max_pos,3)-1))/2,...
       % ID_Node(ID_Element(max_pos,3),3)+Displacement(2*ID_Element(max_pos,3)),...
       % ['最大拉应力',num2str(max_stress)])

%text((ID_Node(ID_Element(min_pos,2),2)+Displacement(2*ID_Element(min_pos,2)-1)+ ...
        %ID_Node(ID_Element(min_pos,3),2)+Displacement(2*ID_Element(min_pos,3)-1))/2,...
       % ID_Node(ID_Element(min_pos,2),3)+Displacement(2*ID_Element(min_pos,2)),...
       % ['最大压应力',num2str(min_stress)])

%%绘制最大竖直位移节点
%subplot(2,2,4);
%title('最大竖直位移图');
%for i=1:Num_Element
  %  axis([-1 110 0 80]);
    %%绘制变形后单元
   % line([ID_Node(ID_Element(i,2),2)+Displacement(2*ID_Element(i,2)-1) ...
       % ID_Node(ID_Element(i,3),2)+Displacement(2*ID_Element(i,3)-1)],...
       % [ID_Node(ID_Element(i,2),3)+Displacement(2*ID_Element(i,2)) ...
       % ID_Node(ID_Element(i,3),3)+Displacement(2*ID_Element(i,3))],'color','r');
   % hold on
%end

%for j=1:Num_Node
     %%绘制变形后节点
   % plot(ID_Node(j,2)+Displacement(2*j-1),ID_Node(j,3)+Displacement(2*j),'k.','MarkerSize',15);
%end

%plot(ID_Node(max_Displacement_pos,2)+Displacement(2*max_Displacement_pos-1), ...
    %ID_Node(max_Displacement_pos,3)+Displacement(2*max_Displacement_pos),'r.','MarkerSize',30);
%text(ID_Node(max_Displacement_pos,2)+Displacement(2*max_Displacement_pos-1), ...
    %ID_Node(max_Displacement_pos,3)+Displacement(2*max_Displacement_pos),['最大竖直位移',num2str(max_Displacement_y)])





