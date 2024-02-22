function plotu(nel,U,nodes,x0,nnode,idisp)
% 显示应力云图
% 输入参数
% iStress 应力分量指标，它可以是下面的值
%           1——u方向位移
%           2——v方向位移  
%确定图形窗口坐标
switch idisp
    case 1
        title='u';
    case 2
        title='v';
end
a=zeros(nnode,2);
for i=1:nel
    a(i,1)=U(2*i-1);
    a(i,2)=U(2*i);
end
%创建图形窗口，并隐藏坐标轴，指定标题
axis equal;
axis off;
set(gcf,'NumberTitle','off');
set(gcf,'Name',title);               
aMin=min(a(:,idisp));  %最小位移
aMax=max(a(:,idisp));  %最大位移
caxis([aMin,aMax]);           %指定颜色映射表
colormap('jet');
%根据单元节点坐标和应力值绘制填充的四边形，显示应力云图
for ie=1:1:nel
    x=[x0(nodes(ie,1),1);
       x0(nodes(ie,2),1);
       x0(nodes(ie,3),1);
       x0(nodes(ie,4),1)];
    y=[x0(nodes(ie,1),2);
       x0(nodes(ie,2),2);
       x0(nodes(ie,3),2);
       x0(nodes(ie,4),2)];
   c=[ a(nodes(ie,1),idisp);
       a(nodes(ie,2),idisp);
       a(nodes(ie,3),idisp);
       a(nodes(ie,4),idisp)];
   set(patch(x,y,c),'EdgeColor','interp')
end
%绘制颜色条，用来指示不同颜色所对应的应力值
yTick=aMin:(aMax-aMin)/10:aMax;
Label=cell(1,length(yTick));
for i=1:length(yTick)
    Label{i}=sprintf('%2f',yTick(i));
end
set(colorbar('vert'),'YTick',yTick,'YTickLabelMode','Manual','YTickLabel',Label);