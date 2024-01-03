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
x_row=find(abs(Node_info(:,2)-2)<1e-9); %找到右上角和右下角的节点标号，在云图中标注数值
y1_row=find(abs(Node_info(:,3)-1)<1e-9);
y2_row=find(abs(Node_info(:,3))<1e-9);
node_num1=intersect(x_row,y1_row);
node_num2=intersect(x_row,y2_row);
text_x2y1=u11(node_num1);
text_x2y2=u11(node_num2);
text_1=['[','(',num2str(Node_info(node_num1,2)),',',num2str(Node_info(node_num1,3)),')',',',num2str(text_x2y1),']'];
text_2=['[','(',num2str(Node_info(node_num2,2)),',',num2str(Node_info(node_num2,3)),')',',',num2str(text_x2y2),']'];
text(Node_info(node_num1,2),Node_info(node_num1,3),text_1,'HorizontalAlignment', 'right');
text(Node_info(node_num2,2),Node_info(node_num2,3),text_2,'HorizontalAlignment', 'right');
subplot(2,1,2);
pcolor(X,Y,u2);
shading interp %色彩平滑
colorbar %分布柱图 
axis equal
title('竖直方向位移分布云图/m')
box on;
x_row=find(abs(Node_info(:,2)-2)<1e-9); %找到右上角和右下角的节点标号，在云图中标注数值
y1_row=find(abs(Node_info(:,3)-1)<1e-9);
y2_row=find(abs(Node_info(:,3))<1e-9);
node_num1=intersect(x_row,y1_row);
node_num2=intersect(x_row,y2_row);
text_x2y1=u22(node_num1);
text_x2y2=u22(node_num2);
text_1=['[','(',num2str(Node_info(node_num1,2)),',',num2str(Node_info(node_num1,3)),')',',',num2str(text_x2y1),']'];
text_2=['[','(',num2str(Node_info(node_num2,2)),',',num2str(Node_info(node_num2,3)),')',',',num2str(text_x2y2),']'];
text(Node_info(node_num1,2),Node_info(node_num1,3),text_1,'HorizontalAlignment', 'right');
text(Node_info(node_num2,2),Node_info(node_num2,3),text_2,'HorizontalAlignment', 'right');
end
