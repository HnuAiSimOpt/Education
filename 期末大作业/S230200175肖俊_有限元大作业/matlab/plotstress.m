function plotstress(nel,stress_node,nodes,x0,iStress)
% 显示应力云图
% 输入参数
% iStress 应力分量指标，它可以是下面的值
%           1——x方向应力
%           2——y方向应力
%           3——剪应力  
%%%%%%%%%%%%%%%%%%%%%%%%%确定图形窗口坐标
switch iStress
    case 1
        title='sigamx';
    case 2
        title='sigamy';
    case 3
        title='sigamxy';
end
%%%%%%%%%%%%%%%%%%%%%创建图形窗口，并隐藏坐标轴，指定标题
axis equal;
axis off;
set(gcf,'NumberTitle','off');
set(gcf,'Name',title);               
stressMin=min(stress_node(iStress,:));  %确定最小应力
stressMax=max(stress_node(iStress,:));  %最大应力
caxis([stressMin,stressMax]);           %指定颜色映射表
colormap('jet');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%根据单元节点坐标和应力值绘制填充的四边形，显示应力云图
for ie=1:1:nel
    x=[x0(nodes(ie,1),1);
       x0(nodes(ie,2),1);
       x0(nodes(ie,3),1);
       x0(nodes(ie,4),1)];
    y=[x0(nodes(ie,1),2);
       x0(nodes(ie,2),2);
       x0(nodes(ie,3),2);
       x0(nodes(ie,4),2)];
   c=[ stress_node(iStress,nodes(ie,1));
       stress_node(iStress,nodes(ie,2));
       stress_node(iStress,nodes(ie,3));
       stress_node(iStress,nodes(ie,4))];
   set(patch(x,y,c),'EdgeColor','interp')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%绘制颜色条，用来指示不同颜色所对应的应力值
yTick=stressMin:(stressMax-stressMin)/10:stressMax;
Label=cell(1,length(yTick));
for i=1:length(yTick)
    Label{i}=sprintf('%2f',yTick(i));
end
set(colorbar('vert'),'YTick',yTick,...
    'YTickLabelMode','Manual','YTickLabel',Label);

