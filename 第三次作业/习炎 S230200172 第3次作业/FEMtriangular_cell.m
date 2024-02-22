% 习炎 S230200172
% 采用 线性三角形单元 实现有限元分析简单 悬臂梁受拉问题
% 板长为0.3m，宽为0.1m，板厚为0.02m。
% 假设与边界条件：左端面固定，右端面受均布载荷w=1000KN/m2，弹性模量为2.1e11Pa，泊松比为0.3。
% 划分6个三角形单元，8个节点，节点编号1，5为固定节点，4，8为受力节点。

%%主程序
clc
clear;

%%初始化
node = double([1   0     0     0
               2   0.1   0     0
               3   0.2   0     0
               4   0.3   0     0
               5   0     0.1   0
               6   0.1   0.1   0
               7   0.2   0.1   0
               8   0.3   0.1   0]);      % 节点信息，第一列为节点编号，2~4列分别为x,y,z方向坐标

ele = double([1   1   2   6
              2   1   5   6
              3   2   3   7
              4   2   6   7
              5   3   4   8
              6   3   7   8]);           % 单元信息，第一列为单元编号，后面各列为单元上的节点号码

num_ele=size(ele, 1);          % 单元数

E = 2.1e11;               % 弹性模量
t = 0.02;                 % 单元厚度
miu = 0.3;                % 泊松比

n_ele = length(ele(:, 1));   %单元数


%组装总体刚度矩阵
dof = length(node(:, 1))*2;       % 自由度数，梁单元每个节点有3个自由度
                                  % (横向位移、扭转角位移、弯曲角位移)
f = ones(dof, 1)*1e8;             % 结构整体外载荷矩阵，整体坐标系下
f_loc = zeros(6, 1);              % 单元外载荷矩阵，局部坐标系下
u = ones(dof, 1)*1e6;             % 位移矩阵
K = zeros(dof);                   % 总体刚度矩阵
stress = zeros(n_ele, 1);         % 单元应力矩阵

for i = 1 : n_ele
    k_ele = TriangleElementStiffness(E, miu, t, node(ele(i, 2:4), 2:4));
    K = assemTriangle(K, k_ele, ele(i, 2), ele(i, 3), ele(i, 4));
end
 
%力边界条件 
f(7)=1.0;        % 4节点横向力
f(8)=0;            % 4节点垂向力
f(15)=1.0;        % 8节点横向力
f(16)=0;            % 8节点垂向力
f(3)=0;f(4)=0;f(5)=0; f(6)=0;
f(11)=0;f(12)=0;f(13)=0;f(14)=0;    

%位移边界条件
u(1)=0; u(2)=0; u(9)=0; u(10)=0; 

%求解未知自由度
index = [];           % 未知自由度的索引
p = [];               % 未知自由度对应的节点力矩阵
for i = 1:dof
    if u(i) ~= 0
        index = [index, i];
        p = [p; f(i)];
    end
end
u(index) = K(index, index) \ p;    % 高斯消去
f = K * u;

% 单元应力
stress = zeros(num_ele, 3);
x1 = node(:, 2) + u(1:2:16);
y1 = node(:, 3) + u(2:2:16);

%可视化
figure;axis equal;
xlim([0,0.4])
ylim([-0.05,0.15])
for i=1 : n_ele               %单元划分
    patch(node(ele(i, 2:4),2), node(ele(i, 2:4),3), ...
        'w', 'FaceColor', 'none', 'LineStyle', '-','EdgeColor', 'r');
end

figure;axis equal;
xlim([0,0.4])
ylim([-0.05,0.15])
for i=1 : n_ele       
    u1 = [u(2*ele(i, 2)-1);
        u(2*ele(i, 2));
        u(2*ele(i, 3)-1);
        u(2*ele(i, 3));
        u(2*ele(i, 4)-1);
        u(2*ele(i, 4))];
    stress(i, :) = TriangleElementStress(E, miu, node(ele(i ,2:4), 2:3), u1, 1)';   % 单元应力计算
    patch(node(ele(i, 2:4), 2), node(ele(i, 2:4), 3), stress(i, 1), 'FaceColor','flat', 'EdgeColor','k');
end
colormap(jet);  % 使用 jet 颜色图
colorbar;  % 显示颜色条


%%计算 元素刚度矩阵，全局刚度矩阵K，应变矩阵函数

function k_ele=TriangleElementStiffness(E,miu,t,node_ele)
    % TriangleElementStiffness 返回具有元素刚度矩阵。
    % 元素刚度矩阵的大小是6 x 6。

    %节点坐标
    x1 = node_ele(1,1);                
    y1 = node_ele(1,2);
    x2 = node_ele(2,1);                
    y2 = node_ele(2,2);
    x3 = node_ele(3,1);                
    y3 = node_ele(3,2);

    A = abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2);   % 单元面积
    a1 = x2*y3 - y2*x3;
    a2 = y1*x3 - x1*y3;
    a3 = x1*y2 - y1*x2;
    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;
    c1 = x3 - x2;
    c2 = x1 - x3;
    c3 = x2 - x1;
    B = 1/2/A * [b1 0 b2 0 b3 0;
                 0 c1 0 c2 0 c3;
                 c1 b1 c2 b2 c3 b3];
    D = E/(1-miu^2) * [1   miu 0;
                       miu 1   0;
                       0   0   (1-miu)/2];     % 平面应力矩阵

    k_ele=t*A*B'*D*B;    % 单元刚度矩阵
    
end

function k_t = assemTriangle(k_t, k_ele, node1 ,node2, node3)
    % assemTriangle 这个函数将平面三角形元素的元素刚度矩阵k组装到全局刚度矩阵K中。
    % 函数在元素刚度矩阵k被组装后返回全局刚度矩阵K。

    d(1:2) = 2*node1 - 1:2*node1;
    d(3:4) = 2*node2 - 1:2*node2;
    d(5:6) = 2*node3 - 1:2*node3;
    for ii = 1 : 6
        for jj = 1 : 6
            k_t(d(ii), d(jj)) = k_t(d(ii), d(jj)) + k_ele(ii, jj);
        end
    end
end

function str=TriangleElementStress(E, miu, node_ele, u1, p)
    % TriangleElementStress 返回元素应力矩阵，平面应力或平面应变选项。
    
    %节点坐标
    x1 = node_ele(1,1);                
    y1 = node_ele(1,2);
    x2 = node_ele(2,1);                
    y2 = node_ele(2,2);
    x3 = node_ele(3,1);                
    y3 = node_ele(3,2);

    A = (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2;   % 单元面积
    a1 = x2*y3 - y2*x3;
    a2 = y1*x3 - x1*y3;
    a3 = x1*y2 - y1*x2;
    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;
    c1 = x3 - x2;
    c2 = x1 - x3;
    c3 = x2 - x1;
    B = 1/2/A * [b1 0  b2 0  b3 0;
                 0  c1 0  c2 0  c3;
                 c1 b1 c2 b2 c3 b3];
    if p == 1
        D = E / (1-miu^2) * [1   miu 0;
                             miu 1   0;
                             0   0   (1-miu)/2];           % 平面应力矩阵
    elseif p == 2
        D = E / (1+miu) / (1-2*miu) * [1-miu miu 0;
                                        miu 1-miu 0;
                                        0 0 (1-2*miu)/2];     % 平面应力矩阵
    end
    str = D*B*u1;   % 单元应力
end
    
    
    
    

