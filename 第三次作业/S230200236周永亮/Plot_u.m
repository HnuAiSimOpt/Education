%% 后处理函数 绘制节点位移云图
function []=Plot_u(u,Node_info)
u1=[];
u2=[];
for i=1:231
 u1=[u1;u(2*i-1)];
 u2=[u2;u(2*i)];
end
u11=u1;
u22=u2; %将 u1、u2 的数值进行保护，后续标注会用到
X1=linspace(min(Node_info(:,2)),max(Node_info(:,2)),100); %将坐标均分成 100 份
Y1=linspace(min(Node_info(:,3)),max(Node_info(:,3)),100);
[X,Y,u1]=griddata(Node_info(:,2),Node_info(:,3),u1,X1',Y1,'v4'); %'v4'代表插值方法为 matlab4 样条函数内插
[X,Y,u2]=griddata(Node_info(:,2),Node_info(:,3),u2,X1',Y1,'v4');
subplot(2,1,1);
pcolor(X,Y,u1);
shading interp %色彩平滑
colorbar %分布柱图 
axis equal
title('水平方向位移分布云图')
box on;
subplot(2,1,2);
pcolor(X,Y,u2);
shading interp %色彩平滑
colorbar %分布柱图 
axis equal
title('竖直方向位移分布云图/m')
box on;
end