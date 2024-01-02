%% 后处理函数 绘制 Mises 应力分布云图
function []=Mise_stress(Node_info,sigma_xx_node,sigma_yy_node,sigma_xxyy_node)
% 先求主应力
sigma_m=[];
for i=1:size(Node_info,1)
 sigma1=(sigma_xx_node(i)+sigma_yy_node(i))/2+0.5*sqrt((sigma_xx_node(i)-sigma_yy_node(i))^2+4*sigma_xxyy_node(i)*sigma_xxyy_node(i));
 sigma2=(sigma_xx_node(i)+sigma_yy_node(i))/2-0.5*sqrt((sigma_xx_node(i)-sigma_yy_node(i))^2+4*sigma_xxyy_node(i)*sigma_xxyy_node(i));
% 求 Von-Mise 等效应力
 sigma_m=[sigma_m;sqrt(sigma1^2+sigma2^2-sigma1*sigma2)];
end
sigma_mm=sigma_m; %将 sigma_m 的数值进行保护，后续标注会用到
X1=linspace(min(Node_info(:,2)),max(Node_info(:,2)),100); %将坐标均分成 100 份
Y1=linspace(min(Node_info(:,3)),max(Node_info(:,3)),100);
[X,Y,sigma_m]=griddata(Node_info(:,2),Node_info(:,3),sigma_m,X1',Y1,'v4'); 
pcolor(X,Y,sigma_m); %绘制云图
colorbar %分布柱图 
axis equal
title('Mises 应力分布云图/Pa')
box on;
end