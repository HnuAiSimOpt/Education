function handle= femesh3D(Element,Node,color)
%绘制三维有限元模型
%Element：单元连接矩阵
%Node：节点坐标
%color: 面表面填充颜色
if nargin < 4
    color = 'c';
end

f = [];
node = [1 2 3 4;
    5 6 7 8
    1 2 6 5;
    2 3 7 6;
    3 7 8 4;
    4 8 5 1];
for i = [1:6]
    f =   [f; Element(:,node(i,:))];
end
f = unique(f,'rows');
handle = patch('Faces',f,'Vertices',Node,'FaceColor',color,'clipping','off');

view([-1 -1 0.3])
set(gca,'ClippingStyle','rectangle','Clipping','off');
axis equal
end











