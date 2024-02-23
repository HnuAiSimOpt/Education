clc;
clear

%导入节点坐标
x0 = load('C:\Users\alein\Desktop\S200200239_FEM\code\xiao_mesh\txt\x0.txt');
%导入abaqus分析得到的节点应力数据
uus = load('C:\Users\alein\Desktop\S200200239_FEM\code\xiao_mesh\txt\uu.txt');
%导入matlab计算得到的节点应力数据
yingli = load('C:\Users\alein\Desktop\S200200239_FEM\code\xiao_mesh\txt\stress_node.txt');
x=x0(:,2);
x=x';
y=x0(:,3);
y=y';
[m n]=size(uus);
u=uus(:,2);
v=uus(:,3);
z=[];
for i=1:m
    z=[z sqrt(u(i)*u(i)+v(i)*v(i))];
end

%%合位移
%插值
[X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x))',linspace(min(y),max(y)),'cubic')
%云图（等高线图），去掉黑色等高线
figure,contourf(X,Y,Z,'LineStyle','none')   
%保持光滑效果
shading flat;
hold on;
colorbar;
title('合位移云图');
%三维曲面
figure,mesh(X,Y,Z)
colorbar;
title('合位移三维云图');

%%x轴应力
%插值
[X,Y,Z1]=griddata(x,y,yingli(:,1),linspace(min(x),max(x))',linspace(min(y),max(y)),'cubic');
%云图（等高线图），去掉黑色等高线
figure,contourf(X,Y,Z1,'LineStyle','none') 
%保持光滑效果
shading flat;
hold on;
colorbar;
title('x轴应力云图');
%三维曲面
figure,mesh(X,Y,Z1)
colorbar;
title('x轴应力三维云图');

%%y轴应力
%插值
[X,Y,Z2]=griddata(x,y,yingli(:,2),linspace(min(x),max(x))',linspace(min(y),max(y)),'cubic');
%云图（等高线图），去掉黑色等高线
figure,contourf(X,Y,Z2,'LineStyle','none')   
%保持光滑效果
shading flat;
hold on;
colorbar;
title('y轴应力云图');
%三维曲面
figure,mesh(X,Y,Z2)
colorbar;
title('y轴应力三维云图');

%%切应力
%插值
[X,Y,Z3]=griddata(x,y,yingli(:,3),linspace(min(x),max(x))',linspace(min(y),max(y)),'cubic');
%云图（等高线图），去掉黑色等高线
figure,contourf(X,Y,Z3,'LineStyle','none')   
%保持光滑效果
shading flat;
hold on;
colorbar;
title('切应力云图');
%三维曲面
figure,mesh(X,Y,Z3)
colorbar;
title('切应力三维云图');