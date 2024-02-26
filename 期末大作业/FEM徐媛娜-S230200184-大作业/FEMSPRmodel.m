% 求解薄钢板的铆接区域在均布压力作用下的空间应力应变问题。
clear all;%清除所有
close all;%关闭所有窗口
clc;
%-----参数设置-----%
% 所有属性采用SI单位.
global L s t E nu p D g;%创建全局变量
L=0.060;% 板半径m
s=0.006;% 连接区域的半径m
t=0.02;%板厚度m
D=7.5e3;%材料的质量密度kg/m3
E=2.06E+011;% 材料杨氏模量.pa
nu = 0.26;%泊松比.
p= 1000;%基本压力.pa
g=9.8;%重力加速度
%-----创建集合偏微分模型------%
% 创建模型.
SPRModel = createpde("structural","static-solid");%创建用于求解空间应力问题的结构模型
gm = multicylinder([s,L],t);%创建带有铆接区域的薄钢板3D计算区域
SPRModel.Geometry = gm;%将创建的gm三维几何图形添加到模型容器中
GeometryModelFigure = ...
    figure('Name','Geometry Model','NumberTitle','off');%绘图，三维模型示意图，标题为Geometry Model
GeometryModelAxes = ...
    axes(GeometryModelFigure,'NextPlot','add ',...%在不清空或重置当前图窗的前提下添加新的图形对象。
         'Box','on',...   %将当前坐标区的 Box 属性设置为 "on" 在坐标区周围显示框轮廓。
         'FontName','Times New Roman','FontSize', 16);%设置字体为新罗马，字号16
GeometryModel = ...%绘制PDE几何模型
    pdegplot(SPRModel,'EdgeLabels','on','FaceLabels','on');%查看边标签和面标签。在图的顶部添加空间以清楚地查看顶部边缘。
xlabel(GeometryModelAxes,'$x/{\rm{(m)}}$','Interpreter' ,'latex');%设置x轴标签，单位，使用 LaTeX 标记解释字符
ylabel(GeometryModelAxes,'$y/{\rm{(m)}}$','Interpreter','latex');%设置y轴标签，单位，使用 LaTeX 标记解释字符
zlabel(GeometryModelAxes,'$z/{\rm{(m)}}$','Interpreter','latex');%设置z轴标签，单位，使用 LaTeX 标记解释字符
title(GeometryModelAxes,'Geometry Model');%设置图片名称Geometry Model

%设置模型的材料参数
structuralProperties(SPRModel,"YoungsModulus",E,...
                              "PoissonsRatio",nu,...
                              "MassDensity",D);%设置几何模型的杨氏模量、泊松比、质量密度
structuralBodyLoad(SPRModel, ...
                  "GravitationalAcceleration",[0;0;g]);%设置模型在z轴方向受重力

%设置边界条件，根据板材实际中所有约束设置边界条件
%structuralBC(SPRModel,"Face",6,"Constraint","fixed");%固定三维模型的面6
structuralBC(SPRModel,"Face",4,"Constraint","fixed");%固定三维模型的面4
structuralBC(SPRModel,"Face",5,"Constraint","fixed");%固定三维模型的面5
structuralBoundaryLoad(SPRModel,...
                      'Face',1,'SurfaceTraction',[0.000,0.000,p]);%指定三维模型面1的牵引力。

% 生成网格
generateMesh(SPRModel,...
             'Hmax',0.400 * s,...%创建一个目标最大元素边缘长度为0.5*s的网格。
             'GeometricOrder','quadratic');%表示二次元的三角形或四面体在其角和边中心具有节点。
MeshFigure = ...
    figure('Name','Mesh','NumberTitle','off');%绘图，图片名称Mesh
MeshAxes = ...
   axes(MeshFigure, 'NextPlot', 'add',...%在不清空或重置当前图窗的前提下添加新的图形对象
                    'Box','on',...
                    'FontName','Times New Roman', 'Fontsize',16);%设置字体为新罗马，字号16
hMesh = ...
    pdeplot3D(SPRModel);%绘图，显示网格
title(MeshAxes, 'Mesh');%设置图片名称Mesh

%求解模型
R = solve(SPRModel);
savefig('Cantilever.fig');%保存

% 计算网格位移幅度.
hFigureDisplacement = ...
    figure('Name','Displacement Magnitude','NumberTitle','off');%绘图网格变形
hAxesDisplacement= ...
    axes(hFigureDisplacement,'NextPlot','add ',...
                             'Box','on',...
                             'FontName','Times New Roman','Fontsize', 16);%设置字体为新罗马，字号16
hDisplacement = ...
    pdeplot3D(SPRModel,"ColorMapData",R.Stress.Magnitude, ...
                       "Deformation",R.Displacement)%绘制网格变形幅度
axis equal;
title(hAxesDisplacement,'Displacement Magnitude');%设置图片名称Displacement Magnitude

%法向应力绘制变形形状
hFigureStress = ...
    figure('Name','Stress-zz','NumberTitle','off');%绘图网格变形
hAxesStress= ...
    axes(hFigureStress,'NextPlot','add ',...
                       'Box','on',...
                       'FontName','Times New Roman','Fontsize', 16);%设置字体为新罗马，字号16
hDisplacement = ...
       pdeplot3D(SPRModel,"ColorMapData",R.Stress.szz, ...%用法向应力的z 分量绘制变形形状
                          "Deformation",R.Displacement)
axis equal;
title(hAxesStress,'Stress-zz');%设置图片名称Stress-zz
                      
%法向应变变形显示
hFigureStrain = ...
    figure('Name','Strain-zz','NumberTitle','off');%绘图网格变形
hAxesStrain= ...
    axes(hFigureStrain,'NextPlot','add ',...
                       'Box','on',...
                       'FontName','Times New Roman','Fontsize', 16);%设置字体为新罗马，字号16
hDisplacement = ...
       pdeplot3D(SPRModel,"ColorMapData",R.Strain.ezz, ...
                          "Deformation",R.Displacement)%用法向应变的z分量绘制变形形状
axis equal;
title(hAxesStrain,'Strain-zz');%设置图片名称Stress-zz                    

%米塞斯应力.
hFigureMises = ...
    figure('Name','von Mises Stress','NumberTitle','off');
hAxesMises = ...
    axes(hFigureMises,'NextPlot','add',...
                      'Box','on',...
                      'FontName','Times New Roman','Fontsize',16);%设置字体为新罗马，字号16
hMises = ...
        pdeplot3D(SPRModel,"ColorMapData",R.VonMisesStress, ...
                           "Deformation",R.Displacement)%用米塞斯应力绘制变形形状
axis equal;
title(hAxesMises,'von Mises stress');%设置图片名称von Mises stress