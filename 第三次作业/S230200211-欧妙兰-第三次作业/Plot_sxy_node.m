%% 后处理函数 绘制剪应力云图
function sigma_xxyy_node=Plot_sxy_node(sigma_xy,Node_info,Ele_info)
sigma_xy_node=[];
for i=1:size(Node_info,1)
 E=[];
 E1=find(abs(Ele_info(:,2)-i)<1e-9);
 E2=find(abs(Ele_info(:,3)-i)<1e-9);
 E3=find(abs(Ele_info(:,4)-i)<1e-9);
 E=[E;E1;E2;E3];
 n=size(E,1);
 sxy_sum=0;
 for j=1:n
 sxy_sum=sxy_sum+sigma_xy(E(j));
 end
 sigma_xy_node=[sigma_xy_node;sxy_sum/n]; 
end
sigma_xxyy_node=sigma_xy_node; %将 sigma_xy_node 的数值进行保护，后续标注会用到
X1=linspace(min(Node_info(:,2)),max(Node_info(:,2)),100); %将坐标均分成 100 份
Y1=linspace(min(Node_info(:,3)),max(Node_info(:,3)),100);
[X,Y,sigma_xy_node]=griddata(Node_info(:,2),Node_info(:,3),sigma_xy_node,X1',Y1,'v4'); 
pcolor(X,Y,sigma_xy_node);
x_row=find(abs(Node_info(:,2)-2)<1e-9); %找到右上角和右下角的节点标号，在云图中标注数值
y1_row=find(abs(Node_info(:,3)-1)<1e-9);
y2_row=find(abs(Node_info(:,3))<1e-9);
node_num1=intersect(x_row,y1_row);
node_num2=intersect(x_row,y2_row);
text_x2y1=sigma_xxyy_node(node_num1);
text_x2y2=sigma_xxyy_node(node_num2);
text_1=['[','(',num2str(Node_info(node_num1,2)),',',num2str(Node_info(node_num1,3)),')',',',num2str(text_x2y1),']'];
text_2=['[','(',num2str(Node_info(node_num2,2)),',',num2str(Node_info(node_num2,3)),')',',',num2str(text_x2y2),']'];
text(Node_info(node_num1,2),Node_info(node_num1,3),text_1,'HorizontalAlignment', 'right');
text(Node_info(node_num2,2),Node_info(node_num2,3),text_2,'HorizontalAlignment', 'right');
shading interp %色彩平滑
colorbar %分布柱图 
axis equal
title('\tau_{xy}分布云图/Pa')
box on;
end