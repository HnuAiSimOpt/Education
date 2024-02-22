% 定义几何参数
L = 1.0;       % 结构长度 (m)
H = 0.5;        % 结构高度 (m)
hole_side = 0.2; % 中间正方形孔的边长 (m)
rect_width = 0.2;  % 右上角和右下角矩形孔的宽度 (m)
rect_height = 0.1; % 右上角和右下角矩形孔的高度 (m)

% 创建主矩形
mainRect = [0, L, L, 0; 0, 0, H, H];

% 创建中间的正方形孔
centerHole = [L/2 - hole_side/2, L/2 + hole_side/2, L/2 + hole_side/2, L/2 - hole_side/2;
              H/2 - hole_side/2, H/2 - hole_side/2, H/2 + hole_side/2, H/2 + hole_side/2];

% 创建右上角的矩形孔
topRightHole = [L - rect_width, L, L, L - rect_width;
                H - rect_height, H - rect_height, H, H];

% 创建右下角的矩形孔
bottomRightHole = [L - rect_width, L, L, L - rect_width;
                   0, 0, rect_height, rect_height];

% 使用patch函数绘制几何形状
figure;
patch('Faces', [1, 2, 3, 4], 'Vertices', mainRect', 'FaceColor', 'blue', 'EdgeColor', 'black'); hold on;
patch('Faces', [1, 2, 3, 4], 'Vertices', centerHole', 'FaceColor', 'white', 'EdgeColor', 'black');
patch('Faces', [1, 2, 3, 4], 'Vertices', topRightHole', 'FaceColor', 'white', 'EdgeColor', 'black');
patch('Faces', [1, 2, 3, 4], 'Vertices', bottomRightHole', 'FaceColor', 'white', 'EdgeColor', 'black');
axis equal;
xlim([0, L]);
ylim([0, H]);
xlabel('Length (m)');
ylabel('Height (m)');
title('2D Geometry with Holes')
