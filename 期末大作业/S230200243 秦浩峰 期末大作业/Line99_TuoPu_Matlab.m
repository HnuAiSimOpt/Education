function top(nelx,nely,volfrac,penal,rmin)
nelx = 100;           % 水平方向上的离散单元数;
nely = 40;              % 竖直方向上的离散单元数;
volfrac = 0.5;        % 容积率,材料体积与设计域体积之比;
penal = 3;              % 惩罚因子，SIMP方法;
rmin = 2;                % 敏度过滤半径，防止出现棋盘格现象;
%初始值
x(1:nely,1:nelx)=volfrac;
loop = 0;
change = 1.;
%当目标函数改变量<=0.01时说明收敛，结束迭代
while change > 0.01
    loop = loop + 1;
    xold = x;
%有限元位移数组[U]
    [U] = FE(nelx,nely,x,penal);
% 计算单元刚度矩阵
    [KE] = lk;
% c是用来储存目标函数的变量    
    c = 0.;
% 遍历设计域矩阵元素
    for ely = 1:nely
        for elx = 1:nelx
            n1 = (nely+1)*(elx-1)+ely;
            n2 = (nely+1)*elx+ely;
            Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
            c = c + x(ely,elx)^penal*Ue'*KE*Ue;
             dc(ely,elx) = -penal*x(ely,elx)^(penal-1)*Ue'*KE*Ue;
        end
    end
 
    %无关网络过滤
    [dc] = check(nelx,nely,rmin,x,dc);
    %采用优化准则(OC)求解当前模型
    [x] = OC(nelx,nely,x,volfrac,dc);
    %过程结果输出
    change = max(max(abs(x-xold)));
    disp(['It.:' sprintf('%4i', loop) 'obj.:' sprintf('%10.4f',c) ...
        'Vol.:' sprintf('%6.3f', sum(sum(x))/(nelx*nely)) ...
        'ch.:' sprintf('%6.3f',change)])
    %过程输出图像
    colormap(gray); imagesc(-x); axis equal; axis tight; axis off; pause(1e-6);
end
   disp(U);%输出位移
%OC优化函数
function [xnew]=OC(nelx,nely,x,volfrac,dc)
        l1 = 0; l2 =100000; move = 0.2;
        while(l2-l1>1e-4)
            lmid = 0.5*(l2+l1);
            xnew = max(0.001,max(x-move,min(1. , min(x+move,x.*sqrt(-dc./lmid)))));
            if sum(sum(xnew))-volfrac*nelx*nely >0
                l1 = lmid;
            else
                l2 = lmid;
            end
        end
        
%无关网络过滤函数
function [dcn] = check(nelx,nely,rmin,x,dc)
dcn = zeros(nely,nelx);
for i = 1:nelx
    for j =1: nely
            sum = 0.0;
        for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
            for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
                    fac = rmin-sqrt((i-k)^2+(j-l)^2);
                    sum = sum+max(0,fac);
                    dcn(j,i) = dcn(j,i)+max(0,fac)*x(l,k)*dc(l,k);
            end
        end
            dcn(j,i)=dcn(j,i)/(x(j,i)*sum);
    end
end

%有限元位移求解程序
function [U]=FE(nelx,nely,x,penal)
[KE]  = lk;
K = sparse(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);

for ely = 1:nely
    for elx=1:nelx
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)*elx+ely;
        edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2;2*n1+1;2*n1+2];
        K(edof,edof)=K(edof,edof)+x(ely,elx)^penal*KE;
    end
end

%载荷施加
F(2*(nelx/2+1)*(nely+1),1)=1;
%固定节点编号
fixeddofs = [2*20+1,2*(20+1) ,2*nelx*(nely+1)+2*(nely/2+1),2*nelx*(nely+1)+2*nely/2+1];
%所有点矩阵编号
alldofs = [1:2*(nely+1)*(nelx+1)];
%求差得自由节点编号
freedofs = setdiff(alldofs,fixeddofs);
%U=F/K矩阵求解
U(freedofs,:     )=K(freedofs,freedofs)\F(freedofs,:  );
U(fixeddofs,:   ) = 0;
%单元刚度矩阵函数
function[KE]=lk
E=1.;
nu=0.3;
k=[1/2-nu/6, 1/8+nu/8, -1/4-nu/12, -1/8+3*nu/8, -1/4+nu/12, -1/8-nu/8, nu/6, 1/8-3*nu/8];
KE = E/(1-nu^2)*[k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8); 
                              k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3);
                              k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2); 
                              k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5);
                              k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4);
                              k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7);
                              k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6);
                              k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1);];

