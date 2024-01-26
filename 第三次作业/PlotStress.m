function PlotStress(iStress)

switch iStress

case 1

title ='x方向正应力';            %图形标题显示“x方向正应力”

case 2

title ='y方向正应力';         %图形标题显示“y方向正应力”

case 3

title ?='剪应力';               %图形标题显示“剪应力”

case 4

title ='最大主应力';           %图形标题显示“最大主应力”

case 5

title ='最小主应力';            %图形标题显示“最小主应力”

 end

figure;                         %创立图形

axis equal ;                 %均分坐标轴

axis off ;                     %关闭坐标轴

set(gcf, 'NumberTitle','off');     %关闭NumberTitle

set(gcf,'Name',title) ;            %图形标题显示title的返回值

 for ie=1:1:16                              %绘制16个单元的相应应力云图

x=[gElementcoordinate(ie,1);gElementcoordinate(ie,3);gElementcoordinate(ie,5)];

y=[gElementcoordinate(ie,2);gElementcoordinate(ie,4);gElementcoordinate(ie,6)];

c=[gElementStress(ie,iStress);gElementStress(ie,iStress);gElementStress(ie,iStress)];

patch(x,y,c)

end

