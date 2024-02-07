% 学号：S230200190，姓名：杜元杰
% 8 节点 6 面体等参单元 求解悬臂梁的应力应变
% 宽度、高度0.2m长度1m，网格数量1080个
% 一端固支，最后一个节点 x 方向受力为2000，y、z方向受力为 1000N，均使用国际单位制
% 位移计算结果与Abqus使用C3D8单元相差很小
format short
clear all;
x_length = 0.2;
y_length = 0.2;
z_length = 1;
lx = 6;%x方向单元个数
ly = lx;%y方向单元个数
lz = 5 * ly;%z方向单元个数
elements_num = lx * ly * lz;%单元个数
nodes_num = (lx + 1) * (ly + 1) * (lz + 1);
elements = zeros(elements_num, 8);%每行为单元包含节点的索引

%计算每个单元包含节点的索引
count = 0;
for i=1:lx
    for j=1:ly
        for k=1:lz
            count = count + 1;
            elements(count,:) =[
             (lx + 1) * (ly + 1) * (k - 1) + (lx +1) * (j - 1) + i,
             (lx + 1) * (ly + 1) * (k - 1) + (lx +1) * (j - 1) + i + 1,
             (lx + 1) * (ly + 1) * (k - 1) + (lx +1) * (j ) + i + 1,
             (lx + 1) * (ly + 1) * (k - 1) + (lx +1) * (j ) + i,
             (lx + 1) * (ly + 1) * (k ) + (lx +1) * (j - 1) + i,
             (lx + 1) * (ly + 1) * (k ) + (lx +1) * (j - 1) + i + 1,
             (lx + 1) * (ly + 1) * (k ) + (lx +1) * (j ) + i + 1,
             (lx + 1) * (ly + 1) * (k ) + (lx +1) * (j ) + i,];
        end
    end
end

%计算每个节点坐标
x_nodes = zeros(nodes_num,1);%节点x坐标
y_nodes = zeros(nodes_num,1);%节点y坐标
z_nodes = zeros(nodes_num,1);%节点z坐标
count = 0;
for k=1:lz+1
    for j=1:ly+1
        for i=1:lx+1
            count = count + 1;
            x_nodes(count) = (i - 1) / lx * x_length;
            y_nodes(count) = (j - 1) / ly * y_length;
            z_nodes(count) = (k - 1) / lz * z_length;
        end
    end
end

%绘图
figure(1)
hold on
axis equal
axis off
plot3(x_nodes,   z_nodes, y_nodes,'o');
xlabel('x');
ylabel('y');
zlabel('z');
view(-45, -30);
for i = 1 : nodes_num
    text(x_nodes(i),   z_nodes(i), y_nodes(i),num2str(i),'FontSize',8);
end
for i = 1 : elements_num
    x = [];
    y = [];
    z = [];
    for j = 1 : 8
        x = [x; x_nodes(elements(i, j))];
        y = [y; y_nodes(elements(i, j))];
        z = [z; z_nodes(elements(i, j))];
    end
    for j = 1 : 4
        line([x(j),x(j + 4)], [z(j),z(j + 4)], [y(j),y(j + 4)] );
    end
    for j = 0 : 4 : 4
        line([x(1+j),x(2+j),x(3+j),x(4+j),x(1+j)] ,[z(1+j),z(2+j),z(3+j),z(4+j),z(1+j)],[y(1+j),y(2+j),y(3+j),y(4+j),y(1+j)]);
    end
end

% 计算三维弹性矩阵D
E = 200e9;% Q235，单位都是国际单位制
poisson = 0.3;% Q235
D = E /((1+poisson)*(1-2*poisson)) *...
    [(1-poisson)  poisson  poisson   0   0    0; 
    poisson  (1-poisson)   poisson   0   0    0;
    poisson  poisson  (1-poisson)    0   0    0;
    0    0    0    (1-2*poisson)/2   0   0;
    0    0    0    0    (1-2*poisson)/2  0;
    0    0    0    0    0   (1-2*poisson)/2];

% 形函数
syms m n o;
N = [1/8 * (1 - m) * (1 - n) * (1 - o);
	 1/8 * (1 + m) * (1 - n) * (1 - o);
	 1/8 * (1 + m) * (1 + n) * (1 - o);
	 1/8 * (1 - m) * (1 + n) * (1 - o);
	 1/8 * (1 - m) * (1 - n) * (1 + o);
	 1/8 * (1 + m) * (1 - n) * (1 + o);
	 1/8 * (1 + m) * (1 + n) * (1 + o);
	 1/8 * (1 - m) * (1 + n) * (1 + o)];
N = N.';

% 8 个积分点等参坐标
xyz = [-1 / sqrt(3), -1 / sqrt(3), -1 / sqrt(3);
       1 / sqrt(3), -1 / sqrt(3), -1 / sqrt(3);
       1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3);
       -1 / sqrt(3), 1 / sqrt(3), -1 / sqrt(3);
       -1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3);
       1 / sqrt(3), -1 / sqrt(3), 1 / sqrt(3);
       1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3);
       -1 / sqrt(3), 1 / sqrt(3), 1 / sqrt(3);];
   
% 计算单元刚度矩阵 k_element 并组装系统刚度矩阵 k_system 
k_system = sparse(3 * nodes_num, 3 * nodes_num);
for i = 1 : elements_num
    x = [];
    y = [];
    z = [];
    index = [];
    for j = 1 : 8
        x = [x; x_nodes(elements(i, j))];
        y = [y; y_nodes(elements(i, j))];
        z = [z; z_nodes(elements(i, j))];
        index = [index;
                 elements(i, j) * 3 - 2;
                 elements(i, j) * 3 - 1;
                 elements(i, j) * 3 ];
    end
    
    % 求3×3 J矩阵
    JJ = zeros(3,3,8);% 8 个点的 J矩阵
    dN = zeros(3,8,8);% 8 个点对等参坐标求导的矩阵
    for j = 1 : 8
        xx = xyz(j, 1);
        yy = xyz(j, 2);
        zz = xyz(j, 3);
        dN(:,:,j) = [double(subs(diff(N,m),{n, o}, {yy, zz}))
        double(subs(diff(N,n),{m, o}, {xx, zz}))
        double(subs(diff(N,o),{m, n}, {xx, yy}))];
        JJ(:, :, j) = dN(:,:,j) * [x, y, z];
    end
        
    %求6×24 B矩阵
    BB = zeros(6,24,8);% 8 个点的B矩阵
    for k = 1:8
        B = [];
        for j = 1 : 8
            dNk = dN(:,:,k);% 第 k 个点的 N 矩阵
            d = JJ(:,:,k) \ dNk(:,j);
            Bi = [d(1) 0 0;% 求 B 的子块
                  0 d(2) 0;
                  0 0 d(3);
                  0 d(3) d(2) ;
                  d(3) 0 d(1);
                  d(2) d(1) 0];
            B = [B,Bi];
        end
        BB(:,:,k) = B;
    end
    
    % 计算单元刚度矩阵 k_element
    k_element = zeros(24, 24);
    
    % 8 个点高斯积分，权重均为 1
    for j = 1 : 8
        k_element = k_element + BB(:, :, j).' * D * BB(:, :, j) * abs(det(JJ(:, :, j)));
    end
    
    % 组装系统刚度矩阵 k_system
    for k = 1 : 24
        for j =  1 : 24
            k_system(index(k), index(j)) = k_system(index(k), index(j)) + k_element(k, j);
        end
    end   
end

% 位移边界条件
bcnode = [];
bcdof = [];
bcval = [];
for i = 1 : (lx + 1) * (ly + 1)
    bcnode = [bcnode; i];%边界条件的节点序号
    bcdof = [bcdof; i * 3 - 2; i * 3 - 1; i * 3];%边界点自由度序号
    bcval = [bcval; 0; 0; 0];%边界点自由度的值
end

% 系统 f 向量，最后一个点x、y、z方向均受力 1000 N
f_system = sparse(nodes_num * 3, 1);
f_system(nodes_num * 3 - 2) = 2000;
f_system(nodes_num * 3 - 1) = 1000;
f_system(nodes_num * 3 ) = 1000;

% 施加边界条件，置一法
for i = 1 : length(bcdof)
    c = bcdof(i);
    for j = 1 : nodes_num * 3
        k_system(c, j) = 0;
    end
    k_system(c, c) = 1;
    f_system(c) = bcval(i);
end

% LU分解 + 高斯消元 计算 位移矩阵，K * d = f
k_system = full(k_system);
[L, U]=lu(k_system);% K = L * U 
utemp = L \ f_system;% utemp = L^-1 * f
d_system = U \ utemp;% d = U^-1 * utemp

% 计算单元应力应变
stress_element = zeros(8, 6, elements_num);% 每个单元内应力不是常数
strain_element = zeros(8, 6, elements_num);% 每个单元内应变不是常数
for i = 1 : elements_num
    x = [];
    y = [];
    z = [];
    index = [];
    for j = 1 : 8
        x = [x; x_nodes(elements(i, j))];
        y = [y; y_nodes(elements(i, j))];
        z = [z; z_nodes(elements(i, j))];
        index = [index;
                 elements(i, j) * 3 - 2;
                 elements(i, j) * 3 - 1;
                 elements(i, j) * 3 ];
    end
    
    %求J矩阵
    JJ = zeros(3,3,8);% 每个点的 J 矩阵
    dN = zeros(3,8,8);% 每个点对等参坐标求导的矩阵
    count = 0;
    for j = 1 : 8
        xx = xyz(j, 1);
        yy = xyz(j, 2);
        zz = xyz(j, 3);
        dN(:,:,j) = [double(subs(diff(N,m),{n, o}, {yy, zz}))
        double(subs(diff(N,n),{m, o}, {xx, zz}))
        double(subs(diff(N,o),{m, n}, {xx, yy}))];
        JJ(:, :, j) = dN(:,:,j) * [x, y, z];
    end
    
    %求B矩阵
    BB = zeros(6,24,8);%每个点的B矩阵
    for k = 1:8
        B = [];
        for j = 1 : 8
            dNk = dN(:,:,k);
            d = JJ(:,:,k) \ dNk(:,j);
            Bi = [d(1) 0 0;%求B的子块
                  0 d(2) 0;
                  0 0 d(3);
                  d(2) d(1) 0;
                  0 d(3) d(2);
                  d(3) 0 d(1)];
            B = [B,Bi];
        end
        BB(:,:,k) = B;
    end
    
   
    d_element = [];% 单元位移矩阵
    %计算8个积分点的应力
    for j = 1 : 24
        d_element = [d_element; d_system(index(j))];
    end
    stress1 = [];% 8个积分点的应力
    for j = 1 : 8
        stress1 = [stress1; (D * BB(:, :, j) * d_element).'];
    end
    
    % 积分点应力->节点应力
    stress2 = zeros(8, 6);% 8个节点的应力
    strain2 = zeros(8, 6);% 8个节点的应变
    for j = 1 : 8;
        xx = xyz(j, 1);
        yy = xyz(j, 2);
        zz = xyz(j, 3);
        NN = double(subs(N, {m, n, o}, {xx, yy, zz}));
        for k = 1 : 8
            % 节点的应力等于 积分点处形函数 * 积分点应力，再求和
            stress2(j, :) = stress2(j, :) + NN(k) * stress1(k, :);
        end
    end
    
    % 节点应变->节点应变
    for j = 1 : 8
        strain2(j, :) = (D \ (stress2(j, :)).').';
    end    
    
    % 保存单元应力和应变
    stress_element(:, :, i) = stress2;
    strain_element(:, :, i) = strain2;
end

% 重合的节点应力和应变取平均值
stress_system = zeros(nodes_num, 6);
strain_system = zeros(nodes_num, 6);
count = zeros(nodes_num,1);% 重复次数
for i = 1 : elements_num
    for j = 1 : 8
        node = elements(i, j);
        stress_system(node, :) = strain_system(node, :) + stress_element(j, :, i);
        strain_system(node, :) = strain_system(node, :) + strain_element(j, :, i);
        count(node, 1) = count(node, 1) + 1;
    end
end
for i = 1 : nodes_num
    stress_system(i,:) = stress_system(i,:) / count(i);
    strain_system(i,:) = strain_system(i,:) / count(i);
end

% 输出数据导techplot
fid_out=fopen('1000.plt','w');
fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
fprintf(fid_out,'VARIABLES="x" "y" "z" "u" "v" "w" "sigax"  "sigmay" "sigmaz" "sigmaxy" "sigmayz" "sigmaxz"\n');
% ET=BRICK，设置为八节点四边形单元
fprintf(fid_out,'ZONE T="flow-field", N= %8d,E=%8d,ET=BRICK, F=FEPOINT\n',nodes_num,elements_num);
d_system = full(d_system);
% 每个节点坐标、位移、应力
for i=1:nodes_num
       fprintf(fid_out,'%16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e \n',x_nodes(i),y_nodes(i), z_nodes(i), d_system(3*i-2), d_system(3*i - 1), d_system(3*i), stress_system(i,1), stress_system(i,2), stress_system(i,3), stress_system(i,4), stress_system(i,5), stress_system(i,6));
end
% 输出每个单元节点的序号
for i=1:elements_num
      fprintf(fid_out,'%8d%8d%8d%8d%8d%8d%8d%8d\n',elements(i,1),elements(i,2),elements(i,3),elements(i,4),elements(i,5),elements(i,6),elements(i,7),elements(i,8));
end 
        

   


