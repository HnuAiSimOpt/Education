function [Node,Element] =generate_mesh(Lx,Ly,Lz,nx,ny,nz)
%网格划分
%输入立方体长宽高（Lx,Ly,Lz）以及每方向对应的网格个数（nx,ny,nz）

p0 = [0,0,0];
x = linspace(0,Lx,nx+1);
y = linspace(0,Ly,ny+1);
z = linspace(0,Lz,nz+1);
temp_X = repmat(x,[length(z),1]);
temp_Y = repmat(y,[length(z),1]);

%调用mesh_lofting函数进行网格划分，得到节点坐标Node和单元连接关系Element
[Node,Element] = mesh_lofting(p0,temp_X,temp_Y,z);

end
