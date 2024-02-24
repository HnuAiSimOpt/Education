% 加载PDE工具箱
pdeModel = createpde();

% 定义主结构和孔的几何形状
R1 = [3,4,0,L,L,0,0,0,H,H]'; % 主矩形
C1 = [3,4,L/2-hole_side/2,L/2+hole_side/2,L/2+hole_side/2,L/2-hole_side/2,H/2-hole_side/2,H/2-hole_side/2,H/2+hole_side/2,H/2+hole_side/2]'; % 中间正方形孔
R2 = [3,4,L-rect_width,L,L,L-rect_width,H-rect_height,H-rect_height,H,H]'; % 右上角矩形孔
R3 = [3,4,L-rect_width,L,L,L-rect_width,0,0,rect_height,rect_height]'; % 右下角矩形孔

% 将几何形状组合成一个集合
gd = [R1,C1,R2,R3];
ns = char('R1','C1','R2','R3')';
sf = 'R1-C1-R2-R3';
dl = decsg(gd,sf,ns);

% 将几何形状添加到PDE模型
geometryFromEdges(pdeModel,dl);

% 生成网格
% mesh = generateMesh(pdeModel,'Hmax',0.03); % 'Hmax' 控制最大单元格尺寸
 mesh = generateMesh(pdeModel, 'Hmax', 0.03, 'GeometricOrder', 'linear');
% 绘制网格
figure;
pdeplot(pdeModel);
axis equal;
title('Mesh of the Geometry with Holes');

