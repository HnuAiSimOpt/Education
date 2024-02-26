function handle = fecontour3D(Element,Node,V)
%绘制三维应力、应变、位移云图
%Axes 为坐标轴
%Element为单元连接矩阵
%Node为节点坐标
%V为节点处的应力/应变/位移值

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

if size(V,1)==1
    V = V';
end

f = unique(f,'rows');
set(gca,'ClippingStyle','rectangle');
handle = patch(gca,'Faces',f,'Vertices',Node,'FaceVertexCData',V,'FaceColor','interp',...
    'clipping','off','LineStyle','-');
view([-1 -1 0.3])
axis equal
colormap(flipud(jet));
colorbar;


set(gca,'ClippingStyle','rectangle','Clipping','off');

end

