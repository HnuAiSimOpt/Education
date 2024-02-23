function plotstrain(nel,stress_node,nodes,x0,iStress)
% 显示应变云图
% 输入参数
% iStress 应力分量指标，它可以是下面的值
%           1——x方向应变
%           2——y方向应变
%           3——剪应变 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%确定图形窗口坐标
format  short 
switch iStress
    case 1
        title='x方向应变';
    case 2
        title='y方向应变';
    case 3
        title='剪应变';
end
%%%%%%%%%%%%%%%%%%%%%创建图形窗口，并隐藏坐标轴，指定标题
figure;
axis equal;
axis off;
set(gcf,'NumberTitle','off');
set(gcf,'Name',title);               
strainMin=min(stress_node(iStress,:));  %确定最小应变
strainMax=max(stress_node(iStress,:));  %最大应变
caxis([strainMin,strainMax]);           %指定颜色映射表
colormap('jet');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%根据单元节点坐标和应力值绘制填充的四边形，显示应变云图
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
yTick=strainMin:(strainMax-strainMin)/10:strainMax;
Label=cell(1,length(yTick));
for i=1:length(yTick)
    Label{i}=sprintf('%2f',yTick(i));
end
set(colorbar('vert'),'YTick',yTick,...
    'YTickLabelMode','Manual','YTickLabel',Label);
