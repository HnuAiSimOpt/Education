%有限元分析矩形薄板，采用四边形单元，一共2000个单元
%姓名：吴靓     学号：S230200201
clc;
clear;
%%
%各种属性定义
length_x = 1;     %矩形长
length_y = 0.4;   %矩形宽
%节点编号及x,y坐标计算
for j = 0:40      %垂直方向41个节点
    for i = 1:51  %水平方向51个节点
        Node.Coord(i+j*51,:) = [(i-1)*length_x/50 j*length_y/40];
    end
end
%网格节点代号编排
for m = 0:39
    for n = 1:50  
        Elem.Def(n+m*50,:) = [n+51*m n+1+51*m n+1+51*(m+1) n+51*(m+1); ];
    end
end
Node.Number = size(Node.Coord,1);  %节点数
[Elem.Number] = size(Elem.Def);    %网格数
Node.fixeddofs = [1 2 102];        %被约束的节点自由度（节点1 x、y,节点51 y)
Node.alldofs = [1:1:2*Node.Number]; %节点所有自由度（x,y)
Node.freedofs = setdiff(Node.alldofs,Node.fixeddofs);  %其他自由节点
Node.Force(:,1) = [2042:2:2142];   %承载节点自由度
Node.Force(:,2) = -100;            %载荷施加
Node.Force(1,2) = -50;
Node.Force(size(Node.Force,1),2) = -50;
Mat.E = 210e9;            %定义弹性模量
Mat.nu = 0.3;             %定义泊松比
Mat.thickness = 0.025;    %定义厚度

%%
%求解刚度矩阵
K = zeros(2*Node.Number);
F = zeros(2*Node.Number,1);
U = zeros(2*Node.Number,1);
for e = 1:Elem.Number
    [x] = Node.Coord(Elem.Def(e,:),:);
    X = reshape([x]',1,8);
    syms s t;
    a = (X(2)*(s-1)+X(4)*(-1-s)+X(6)*(1+s)+X(8)*(1-s))/4;
    b = (X(2)*(t-1)+X(4)*(1-t)+X(6)*(1+t)+X(8)*(-1-t))/4;
    c = (X(1)*(t-1)+X(3)*(1-t)+X(5)*(1+t)+X(7)*(-1-t))/4;
    d = (X(1)*(s-1)+X(3)*(-1-s)+X(5)*(1+s)+X(7)*(1-s))/4;
    B1 = [a*(t-1)/4-b*(s-1)/4 0 ; 0 c*(s-1)/4-d*(t-1)/4 ; c*(s-1)/4-d*(t-1)/4 a*(t-1)/4-b*(s-1)/4];
    B2 = [a*(1-t)/4-b*(-1-s)/4 0 ; 0 c*(-1-s)/4-d*(1-t)/4 ; c*(-1-s)/4-d*(1-t)/4 a*(1-t)/4-b*(-1-s)/4];
    B3 = [a*(t+1)/4-b*(s+1)/4 0 ; 0 c*(s+1)/4-d*(t+1)/4 ; c*(s+1)/4-d*(t+1)/4 a*(t+1)/4-b*(s+1)/4];
    B4 = [a*(-1-t)/4-b*(1-s)/4 0 ; 0 c*(1-s)/4-d*(-1-t)/4 ; c*(1-s)/4-d*(-1-t)/4 a*(-1-t)/4-b*(1-s)/4];
    BM = [B1 B2 B3 B4];
    Jdet1 = [0 1-t t-s s-1 ; t-1 0 s+1 -s-t ; s-t -s-1 0 t+1 ; 1-s s+t -t-1 0];
    J = [X(1) X(3) X(5) X(7)]*Jdet1*[X(2);X(4);X(6);X(8)]/8;
    B = BM/J;
    D = (Mat.E/(1-Mat.nu*Mat.nu))*[1,Mat.nu,0;Mat.nu,1,0;0,0,(1-Mat.nu)/2];
    BD = J*B'*D*B;
    r = int(int(BD,t,-1,1),s,-1,1);
    z = Mat.thickness*r;
    KE = double(z);
    node1 = Elem.Def(e,1); 
    node2 = Elem.Def(e,2); 
    node3 = Elem.Def(e,3); 
    node4 = Elem.Def(e,4);
    fdof = [2*node1-1;2*node1;2*node2-1;2*node2;2*node3-1;2*node3;2*node4-1;2*node4];
    K(fdof,fdof) = K(fdof,fdof)+KE;
end

%%
%求解位移及应力
F(Node.Force(:,1)) = Node.Force(:,2);
U(Node.freedofs,:) = K(Node.freedofs,Node.freedofs) \ F(Node.freedofs,:);     
U(Node.fixeddofs,:) = 0;
F=K*U;
for e=1:Elem.Number
    node1 = Elem.Def(e,1); 
    node2 = Elem.Def(e,2); 
    node3 = Elem.Def(e,3); 
    node4 = Elem.Def(e,4);
    fdof = [2*node1-1;2*node1;2*node2-1;2*node2;2*node3-1;2*node3;2*node4-1;2*node4];
    u = U(fdof(:,1),1);
    for i=1:4
        stress = D*B*u;
        L=[-1 -1;1 -1;1 1;-1 -1];
        Jointsigma = subs(stress, {s,t}, L(i,:));
        sigma(fdof([2*i-1 2*i]))= double(Jointsigma([1 2],:));
    end
end

%%
% 绘制位移以及应力云图
x=Node.Coord(1:51)';
y=[0:0.01:0.4];

for i=0:40
    xsigma(i+1,:)=sigma(102*i+1:2:102*i+101);
    ysigma(i+1,:)=sigma(102*i+2:2:102*i+102);
    Uxdisplay(i+1,:)=U(102*i+1:2:102*i+101);
    Uydisplay(i+1,:)=U(102*i+2:2:102*i+102);
end

figure(1)
contourf(x,y,ysigma)
set(gcf,'position',[200,200,800,800])
colormap(hsv);
colorbar;
axis equal;
grid on;
set(gca, 'GridAlpha', 0.5); 
set(gca, 'XTick', 0:1/50:1); 
set(gca, 'YTick', 0:0.4/40:0.4);
ylabel('y');
xlabel('x');
title('stress-y distribution');

figure(2)
contourf(x,y,xsigma)
set(gcf,'position',[200,200,800,800])
colormap(hsv);
colorbar;
axis equal;
grid on;
set(gca, 'GridAlpha', 1); 
set(gca, 'XTick', 0:1/50:1); 
set(gca, 'YTick', 0:0.4/40:0.4);
ylabel('y');
xlabel('x');
title('stress-x distribution');

figure(3)
contourf(x,y,Uydisplay)
set(gcf,'position',[200,200,800,800])
colormap(hsv);
colorbar;
axis equal;
grid on;
set(gca, 'GridAlpha', 1); 
set(gca, 'XTick', 0:1/50:1); 
set(gca, 'YTick', 0:0.4/40:0.4);
ylabel('y');
xlabel('x');
title('Uy distribution');

figure(4)
contourf(x,y,Uxdisplay)
set(gcf,'position',[200,200,800,800])
colormap(hsv);
colorbar;
axis equal;
grid on;
set(gca, 'GridAlpha', 1); 
set(gca, 'XTick', 0:1/50:1); 
set(gca, 'YTick', 0:0.4/40:0.4);
ylabel('y');
xlabel('x');
title('Ux distribution');