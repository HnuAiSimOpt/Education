%学号：S230200190，姓名：杜元杰
%平面三角形单元程序
%6、8节点固定，1节点收到向下的力，大小为1，均使用国际单位制
clear all;
clc;
elements_num = 10;%单元个数
elements = [];%每行为单元包含节点的序号
x_length = 5;
y_length = 1;
it = elements_num / 2;%循环次数
for i=1:2:(it*2)-1
    a = [i, i+1, i+3];
    b = [i, i+2, i+3];
    elements = [elements; a; b];
end
x_nodes = [];%节点x坐标
y_nodes = [];%节点y坐标
for i=1:it + 1
    x_nodes = [x_nodes; (i - 1) / it * x_length; (i - 1) / it * x_length]; 
    y_nodes = [y_nodes; 0; -y_length];
end

%绘图
figure(1)
hold on
axis off %隐藏坐标系
axis equal %坐标轴长度等
for i = 1 : elements_num
         x1 = x_nodes(elements(i, 1));
         x2 = x_nodes(elements(i, 2));
         x3 = x_nodes(elements(i, 3));
         y1 = y_nodes(elements(i, 1));
         y2 = y_nodes(elements(i, 2));
         y3 = y_nodes(elements(i, 3));
         line([x1,x2],[y1,y2]);
         line([x1,x3],[y1,y3]);
         line([x2,x3],[y2,y3]);
end
for i = 1 : 12
    text(x_nodes(i), y_nodes(i) + 0.1, num2str(i));
end

%6,8节点固定
bcnode = [6; 8];%边界条件的节点序号
bcdof = [];%边界点自由度序号
bcval = [];%边界点自由度的值
for i = 1: length(bcnode)
index = bcnode(i);
bcdof = [bcdof; 1+2*(index-1); 2+2*(index-1)];
bcval = [bcval; 0; 0];
end
systemdof = length(x_nodes) * 2;%系统自由度
elementdof = 2 * 3;%每个单元自由度
f_system = sparse(systemdof,1);%系统载荷向量
k = sparse(elementdof,elementdof);%每个单元的k矩阵		
k_system = sparse(systemdof,systemdof);%系统的k矩阵		
d_system = sparse(systemdof,1);%系统的d矩阵		
element_d = sparse(elementdof,1);%每个单元的d矩阵		
stress = zeros(elements_num,3);%应力矩阵		
strain = zeros(elements_num,3);%应变矩阵
index = sparse(elementdof,1);%单元刚度矩阵自由度序号索引
B = sparse(3,elementdof);% 应变矩阵节点位移->单元应变		
D = sparse(3,3);%本构矩阵/弹性矩阵 应变->应力
p = 0;%泊松比
E = 1;%杨氏模量
D = [ 1 p 0;
      p 1 0;
      0 0 (1 - p) / 2];
D = D * E / (1 - p^2);%弹性矩阵
h = 1;

%计算每个单元的刚度矩阵并组装
for i = 1 : elements_num 
    node1 = elements(i, 1);
    node2 = elements(i, 2);
    node3 = elements(i, 3);
    x1 = x_nodes(node1);
    x2 = x_nodes(node2);
    x3 = x_nodes(node3);
    y1 = y_nodes(node1);
    y2 = y_nodes(node2);
    y3 = y_nodes(node3);
    A = 0.5 *((x2 * y3 - x3 * y2) + (y2 - y3) * x1 + (x3 - x2) * y1);
    B = [y2-y3 0 y3-y1 0 y1-y2 0;
         0 x3-x2 0 x1-x3 0 x2-x1;
         x3-x2 y2-y3 x1-x3 y3-y1 x2-x1 y1-y2];
    B = 1 / (2 * A) * B;
    k = h * A * transpose(B) * D * B;%计算单元刚度矩阵
    
    %计算每一行/列的自由度总体序号
    index = [node1 * 2 - 1;
         node1 * 2;
         node2 * 2 - 1;
         node2 * 2 ;
         node3 * 2 - 1;
         node3 * 2;];
    %遍历单元刚度矩阵，加入总体刚度矩阵 
    for m = 1 : 6
        for n = 1 : 6
            k_system(index(m), index(n)) = k_system(index(m), index(n)) + k(m, n);
        end
    end 
end
det_k = det(k_system);
%施加边界条件，置一法
%不施加边界条件 -> k 矩阵奇异
for i = 1 : length(bcdof)
    c = bcdof(i);
    for j = 1 : systemdof
        k_system(c, j) = 0;
    end
    k_system(c, c) = 1;
    f_system(c) = bcval(i);
end

%f向量，节点1向下的力，大小为1
f_system(1) = -1;

%位移矩阵，K d = f
[L U]=lu(k_system);% K = L U 
utemp = L \ f_system;% y = L^-1 * f
det = det(U);

d_system = U \ utemp;% d = U^-1 * y
L = full(L);
U = full(U);
rcond(U)
k_system = full(k_system);
%计算每个单元的节点应力
for i = 1 : elements_num 
    node1 = elements(i, 1);
    node2 = elements(i, 2);
    node3 = elements(i, 3);
    x1 = x_nodes(node1);
    x2 = x_nodes(node2);
    x3 = x_nodes(node3);
    y1 = y_nodes(node1);
    y2 = y_nodes(node2);
    y3 = y_nodes(node3);
    A = 0.5 *((x2 * y3 - x3 * y2) + (y2 - y3) * x1 + (x3 - x2) * y1);
    B = [y2-y3 0 y3-y1 0 y1-y2 0;
         0 x3-x2 0 x1-x3 0 x2-x1;
         x3-x2 y2-y3 x1-x3 y3-y1 x2-x1 y1-y2];
    B = 1 / (2 * A) * B;
    
    %计算每一行/列的自由度总体序号
    index = [node1 * 2 - 1;
         node1 * 2;
         node2 * 2 - 1;
         node2 * 2 ;
         node3 * 2 - 1;
         node3 * 2;];
    %遍历单元刚度矩阵，加入总体刚度矩阵 
    for j = 1 : 6
        element_d(j) = d_system(index(j));
    end
    strain_element = B * element_d;% 计算单元应变
    stress_element = D * strain_element;% 计算单元应力
    for k=1:3
        %三角形单元为常应变单元，应力结果在一个单元内也为常数。
        strain(i,k) = strain_element(k);% 总体单元的应变
        stress(i,k) = stress_element(k); % 总体单元的应力         
    end 
end

%重合的节点应力和应变取平均值
stress_averaged = zeros(12,3);
strain_averaged = zeros(12,3);
count = zeros(12,1);%重复次数
%遍历每个单元，取出单元内节点应力、应变
for i = 1 : elements_num
    for j = 1 : 3
        node = elements(i , j);
        stress_averaged(node, :) = stress_averaged(node, :) + stress(i, :);
        strain_averaged(node, :) = strain_averaged(node, :) + strain(i, :);
        count(node, 1) = count(node, 1) + 1;
    end
end
%取均值
for i = 1 : 12
    stress_averaged(i,:) = stress_averaged(i,:) / count(i);
    strain_averaged(i,:) = strain_averaged(i,:) / count(i);
end
for i = 1 : 12
    fprintf("第 %d 个节点应力：\n σxx = %f，σyy = %f, σxy = %f\n", i, stress_averaged(i, 1), stress_averaged(i, 2), stress_averaged(i, 3));
end

%输出数据
fid_out=fopen('beam01.plt','w');

fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
fprintf(fid_out,'VARIABLES="x" "y" "u" "v" "sigmax"  "sigmay" "sigmaxy"\n');
% ET=TRIANGLE，三角形单元
fprintf(fid_out,'ZONE T="flow-field", N= %8d,E=%8d,ET=TRIANGLE, F=FEPOINT\n',12,10);
d_system = full(d_system);
for i=1:12
       fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',x_nodes(i),y_nodes(i), d_system(2*i-1),d_system(2*i),stress_averaged(i,1),stress_averaged(i,2),stress_averaged(i,3));
end
for i=1:10
      fprintf(fid_out,'%8d%8d%8d\n',elements(i,1),elements(i,2),elements(i,3));
end

    