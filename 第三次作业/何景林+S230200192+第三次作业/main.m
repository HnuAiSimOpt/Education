%% comments：
%{
     姓名：何景林     学号：S230200192
purpose：
考虑一个悬臂梯形梁的平面应力问题，假设梯形下底a=20mm，上底b=10mm,梯形高h=60mm，厚度t=5mm，弹性模量为200000MPa，泊松比为0.3，
用有限元法(平面三角形单元)划分250个单元、156个节点，节点156处施加集中荷载P为1000N，梯形的下底，即编号为1、2、3、4、5、6的三个节点为固定端,求出节点位移，各单元应力和应变。

计算结果见word
%}
%% 1.主函数
%% 1.1输入结点数、单元数、结点自由度数；节点坐标、单元中结点间关系、单元材料和几何参数、边界条件、节点载荷
clear;clc
%% 节点&单元数
nnd = 156 ;               	% 节点数 				
nel = 250 ;                	% 单元数
nelx = 25;                  % 横向均分数
nely = 5;                   % 纵向均分数
lengthx = 60;               % 梯形下底边
lengthy = 20;               % 梯形高
nne = 3 ;                	% 单元节点数			
nodof =2;                	% 节点自由度
eldof = nne*nodof;       	% 单元自由度
%% 单元材料和几何参数
E = 200000.;     % E（Mpa）
vu = 0.3;        % 泊松比
thick = 5.;      % 横梁厚度（mm） 
%% 节点坐标
% geom(nnd,2)存储节点的x,y坐标
geom = [];
for i=1:nelx+1
    for j=1:nely+1
          geom=[geom; (i-1)*lengthx/nelx      -0.5*lengthy*(1+(nelx+1-i)/nelx)*(1-(j-1)/nely)];
    end
end

%% 节点关联关系
%connec(nel,2)存储节点关联关系。【节点关联关系矩阵】
connec=[];
for i=1:nelx
    for j=1:nely
        connec=[connec; (nely+1)*(i-1)+j (nely+1)*i+j (nely+1)*(i-1)+j+1;];
        connec=[connec; (nely+1)*i+j (nely+1)*i+j+1 (nely+1)*(i-1)+j+1;];
    end
end
%% 形成平面应力的弹性矩阵
dee = formdsig(E,vu); 
%% 边界条件
%nf(nnd,nodof)储存边界条件。0代表自由度被约束，1（编号）代表没有。
%%1.假设每个节点每个自由度都未被约束。
nf = ones(nnd, nodof);    % 将矩阵 nf 初始化为 1
%%2.对实际有约束的自由度设置为0
nf(1,1) = 0; nf(1,2) = 0;  % 节点1、2、3、4、5、6为固端约束，自由度0
nf(2,1) = 0; nf(2,2) = 0;  
nf(3,1) = 0; nf(3,2) = 0;  
nf(4,1) = 0; nf(4,2) = 0;  
nf(5,1) = 0; nf(5,2) = 0;  
nf(6,1) = 0; nf(6,2) = 0;
%%3.对未被约束的节点进行编号
n=0;
for i=1:nnd
    for j=1:nodof
        if nf(i,j) ~= 0 
            n=n+1;
           nf(i,j)=n;
        end
    end
end
 
%% 节点载荷
%load矩阵储存荷载
Nodal_loads= zeros(nnd, 2);
Nodal_loads(156,1) = 0.; Nodal_loads(156,2) = -1000.;    % Node 2
 
%% 1.2整体力矢和总刚度组装 
% 整体力矢的组装
fg=zeros(n,1);
for i=1: nnd
    if nf(i,1) ~= 0
       fg(nf(i,1))= Nodal_loads(i,1);
    end
    if nf(i,2) ~= 0
       fg(nf(i,2))= Nodal_loads(i,2);
    end
end
% 总刚的组装
kk = zeros(n, n);
for i=1:nel
    [bee,g,A] = elem_T3(i,nne,nodof,geom,connec,nf);  % 形成应变矩阵bee，操作矢量g，单元面积A
    ke=thick*A*bee'*dee*bee; % 计算单刚
    [kk]=form_KK(kk,ke, g,eldof);    % 组装总刚
end
 
%% 1.3分析
%1 整体有限元方程的求解
delta = kk\fg ; % 求解未知位移
node_disp=zeros(nnd,2);
for i=1: nnd %
    if nf(i,1) == 0 
        x_disp =0.; 
    else
        x_disp = delta(nf(i,1)); 
    end
    if nf(i,2) == 0 
        y_disp = 0.; 
    else
        y_disp = delta(nf(i,2)); 
    end
    node_disp(i,:) =[x_disp y_disp];
end
%2 单元应力和应变
for i=1:nel
    [bee,g,A] = elem_T3(i,nne,nodof,geom,connec,nf); % 形成应变矩阵bee，操作矢量g，单元面积A
    eld=zeros(eldof,1);     % 将单元位移初始化为零
    for m=1:eldof
        if g(m)==0
            eld(m)=0;
        else
            eld(m)=delta(g(m)); % 检索单元位移
        end
    end
    eps=bee*eld;   % 计算应变（应变矩阵*位移）
    EPS(i,:)=eps ; % 存储所有结点的应变
    sigma=dee*eps; % 计算应力
    SIGMA(i,:)=sigma ; % 存储所有节点的应力
end
 
%% 1.4 输出计算结果
fid =fopen('CST_COARSE_MESH_RESULTS.txt','w');%命名计算结果文件
fprintf(fid, '-------------------------------------------------------- \n');
fprintf(fid, ' \n ********** 分析结果 **********\n');
 
%% 打印节点位移
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, '                    节点位移 \n');
fprintf(fid, 'Node      disp_x          disp_y \n');
for i=1:nnd
    fprintf(fid,' %g,     %8.5e,     %8.5e\n',  ...
        i, node_disp(i,1), node_disp(i,2));
end
fprintf(fid,'\n');
 
%% 打印节点应力
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, '                    单元应力 \n');
fprintf(fid, 'element    sigma_(xx)         sigma_(yy)         tau_(xy)\n');
for i=1:nel
    fprintf(fid,' %g,      %7.4e,       %7.4e,       %7.4e\n',i, ...
                SIGMA(i,1),SIGMA(i,2),SIGMA(i,3));
end
 
%% 打印节点应变
fprintf(fid, '------------------------------------------------------ \n');
fprintf(fid, '                     单元应变 \n');
fprintf(fid, 'element    epsilon_(xx)       epsilon_(yy)         gamma_(xy)\n');
for i=1:nel
    fprintf(fid,' %g,      %7.4e,       %7.4e,       %7.4e\n',i, ...
                EPS(i,1),EPS(i,2),EPS(i,3));
end
type CST_COARSE_MESH_RESULTS.txt

% 绘制 x_direction 上的应力（梁的纵向应力云图）
x_stress = SIGMA(:,1);
cmin = min(x_stress);
cmax = max(x_stress);
caxis([cmin cmax])    %设置颜色范围
patch('Faces', connec, 'Vertices', geom, 'FaceVertexCData',x_stress, 'Facecolor','flat','Marker','o')
%patch('Faces',F,'Vertices',V，'FaceVertexCData',x_stress) 
% 创建一个或多个多边形，其中 V 指定顶点的值，F 定义要连接的顶点
% 'FaceVertexCData'这个东西是来指定顶点vertex上的颜色的。
% 'Facecolor','flat'设置属性
title('\sigma_{xx}(KN/mm^{2})','fontname','Times New Roman')
colorbar

%% 2.子函数
%% 2.1单元分析
function[bee,g,A] = elem_T3(i,nne,nodof,geom,connec,nf)
% 该函数返回单元 i 的节点坐标、单元面积、应变矩阵、操作向量
% 返回节点坐标
x1 = geom(connec(i,1),1);   y1 = geom(connec(i,1),2);
x2 = geom(connec(i,2),1);   y2 = geom(connec(i,2),2);
x3 = geom(connec(i,3),1);   y3 = geom(connec(i,3),2);
% 单元面积的计算
A = (0.5)*det([1   x1    y1; 
               1   x2    y2; 
               1   x3    y3]);
 
 m11 = (x2*y3 - x3*y2)/(2*A);
 m21 = (x3*y1 - x1*y3)/(2*A);
 m31 = (x1*y2 - y1*x2)/(2*A);
 m12 = (y2 - y3)/(2*A);
 m22 = (y3 - y1)/(2*A);
 m32 = (y1 - y2)/(2*A);
 m13 = (x3 - x2)/(2*A);
 m23 = (x1 - x3)/(2*A);
 m33 = (x2 -x1)/(2*A);
 
% 应变矩阵的计算
bee = [ m12   0    m22   0    m32      0; ...
         0   m13    0   m23    0     m33; ...
        m13  m12   m23  m22   m33    m32] ;
 
% 操作向量
l=0;
for k=1:nne   
    for j=1:nodof
    l=l+1;
    g(l)=nf(connec(i,k),j);
    end
end
end

%% 2.2.总刚度矩阵组装
function[kk]=form_KK(kk,ke, g,eldof)
% 此函数组装总刚 
for i=1:eldof
   if g(i) ~= 0
      for j=1: eldof
          if g(j) ~= 0
             kk(g(i),g(j))= kk(g(i),g(j)) + ke(i,j);
          end
       end
   end
end
end

%% 2.3.弹性矩阵构建
function[dee] = formdsig(E,vu)
% 此函数形成平面应力问题的弹性矩阵
c=E/(1.-vu*vu);
dee=c*[1          vu         0.        ;...
     vu         1           0.        ;...
     0.          0.       .5*(1.-vu)];
end