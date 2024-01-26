%平面问题常应变三角形单元
%用线性三角形单元计算悬臂梁在一定集中载荷下的节点位移
%雷锐 S230200246
global nnd nel nne nodof eldof n geom connec B D U nf loads x1 y1 x2 y2 x3 y3 ae   %定义全局变量
format short e
%确定节点数和单元数
nnd=24;
nel=30;
nne=3;
nodof=2;
eldof=nne*nodof;
%确定节点坐标矩阵
geom=[0,0;0,10;0,20;0,30;  %第一列节点坐标
    10,0;10,10;10,20;10,30;  %第二列
    20,0;20,10;20,20;20,30;  %第三列
    30,0;30,10;30,20;30,30;  %第四列
    40,0;40,10;40,20;40,30;  %第五列
    50,0;50,10;50,20;50,30;]  %第六列
%确定节点之间的连接关系，逆时针编号
connec=[1,5,2;
    5,6,2;
    2,6,3;
    6,7,3;
    3,7,4;
    7,8,4;
    5,9,6;
    9,10,6;
    6,10,7;
    10,11,7;
    7,11,8;
    11,12,8;
    9,13,10;
    13,14,10;
    10,14,11;
    14,15,11;
    11,15,12;
    15,16,12;
    13,17,14;
    17,18,14;
    14,18,15;
    18,19,15;
    15,19,16;
    19,20,16;
    17,21,18;
    21,22,18;
    18,22,19;
    22,23,19;
    19,23,20;
    23,24,20]   %30个三角形单元的节点连接信息
%确定载荷矩阵
loads=zeros(nnd,2);
loads(24,2)=-1000;   %在最右上角的24号节点有一个竖直向下的大小为1KN的集中载荷
%材料的几何和物理参数
E=2e5;   %单位MPa
v=0.3;   %泊松比0.3
d=5;   %横梁厚度为5mm
%根据公式算出平面应力的弹性矩阵
D=E/(1-v*v)*[1,v,0;v,1,0;0,0,(1-v)/2];
%设置边界条件，定义每个节点的自由度
%0表示该节点该方向上的自由度为0，被约束了；1表示自由度为1，未被约束
nf=ones(nnd,nodof);
for i=1:4
    for j=1:2
        nf(i,j)=0;
    end
end
%靠墙的四个节点被约束，自由度为0
%对其他没有被约束的节点按顺序进行编号
n=0;
for i=5:nnd
    for j=1:nodof
        n=n+1;
        nf(i,j)=n;
    end
end
%整体力矢的组装
fg=zeros(n,1);
for i=5:nnd
        fg(nf(i,1))=loads(i,1);
        fg(nf(i,2))=loads(i,2);
end
%计算并组装总体刚度矩阵
kk=zeros(n,n);
for i=1:nel
    %计算单元节点坐标
    x1=geom(connec(i,1),1); y1=geom(connec(i,1),2);
    x2=geom(connec(i,2),1); y2=geom(connec(i,2),2);
    x3=geom(connec(i,3),1); y3=geom(connec(i,3),2);
    %通过三个节点的坐标计算三角形单元面积
    ae=[1,x1,y1;1,x2,y2;1,x3,y3];
    A=(0.5)*det(ae);
    b1=(y2-y3)/(2*A);
    b2=(y3-y1)/(2*A);
    b3=(y1-y2)/(2*A);
    c1=(x3-x2)/(2*A);
    c2=(x1-x3)/(2*A);
    c3=(x2-x1)/(2*A);
    %组装成应变矩阵
    B=[b1 0 b2 0 b3 0;
        0 c1 0 c2 0 c3;
        c1 b1 c2 b2 c3 b3];
    %操作向量
    l=0;
    for k=1:nne
        for j=1:nodof
            l=l+1;
            g(l)=nf(connec(i,k),j);
        end
    end
    ke=d*A*B'*D*B;  %计算单元刚度矩阵
    %组装系统刚度矩阵
    for i=1:eldof
        if g(i)~=0
            for j=1:eldof
                if g(j)~=0
                    kk(g(i),g(j))=kk(g(i),g(j))+ke(i,j);
                end
            end
        end
    end
end
    %求解矩阵方程
    U=kk\fg;
    for i=1:nnd
        if nf(i,1)==0
            x_disp=0;
        else
            x_disp=U(nf(i,1));
        end
        if nf(i,2)==0
            y_disp=0;
        else
            y_disp=U(nf(i,2));
        end
        node_disp(i,:)=[x_disp y_disp];
    end
    %输出结果
    fprintf('------------------------');
    fprintf('\n****** 输出节点位移如下 ******\n');
    fprintf('------------------------\n');
    %输出节点位移
    fprintf('Node      disp_x          disp_y \n');
    for i=1:nnd
        fprintf(' %g,     %8.5e,     %8.5e\n',  ...
            i, node_disp(i,1), node_disp(i,2));
    end
    
    