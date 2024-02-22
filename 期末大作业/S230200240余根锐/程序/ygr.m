%-------------------------------------------------------------------------
% 三维楔形槽有限元分析 
% 六面体单元
% 余根锐   机械工程2301班  S230200240
%--------------------------------------------------------------------------
% 主要变量描述：
%  x0 = 节点坐标值
% nodes = 单元的节点编号
% nel    % 元素总量 
% KE= 单元刚度矩阵
% K = 系统刚度矩阵
%  F= 系统载荷向量
% disp = 系统节点位移向量
% stress = 应力矩阵
% strain = 应变矩阵
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
%前处理
% 建立模型，网格划分
%-------------------------------------------------------------------------
function main()
ly=8;             % X轴的元素总量
lx=20;            %Y轴的元素总量
lz=8;             %Z轴的元素总量
lengthx=160;
lengthy=40;
lengthz=40;    %单元尺寸
nel=lx*ly*lz;            % 元素总量 
t=10;
nnel=8;
K=sparse(3*(ly+1)*(lx+1)*(lz+1),3*(ly+1)*(lx+1)*(lz+1));%定义系统刚度矩阵 
F=sparse(3*(ly+1)*(lx+1)*(lz+1),1);%定义力矩阵 
disp=zeros(3*(ly+1)*(lx+1)*(lz+1),1); %定义位移矩阵 l=8; 每个元素的节点数   
%--------------------------
% 定义材料
E=2100000; %弹性模量E   
NU=0.3;    %泊松比μ
%--------------------------

%--------------------------
%%节点坐标
%--------------------------
grad=-0.4; 
x0=[];
for k=1:lz+1
   for i=1:lx+1
      if i<=lx/2+1 
            for j=1:ly+1
            x0=[x0; (i-1)*lengthx/lx  (j-1)*(lengthy+(i-1)*lengthx*grad/lx)/ly (k-1)*lengthz/lz   ]; 
            end
      else
            for j=1:ly+1
               x0=[x0; (i-1)*lengthx/lx  (j-1)*(lengthy+0.5*lengthx*grad-(i-lx/2-1)*lengthx*grad/lx)/ly  (k-1)*lengthz/lz ]; 
            end
      end
    end
end
%%节点编号
nodes=[];
for k=1:lz
   for i=1:lx
      for j=1:ly
       nodes=[nodes; (k-1)*(lx+1)*(ly+1)+(ly+1)*(i-1)+j (k-1)*(lx+1)*(ly+1)+(ly+1)*i+j (k-1)*(lx+1)*(ly+1)+(ly+1)*i+j+1 (k-1)*(lx+1)*(ly+1)+(ly+1)*(i-1)+j+1  k*(lx+1)*(ly+1)+(ly+1)*(i-1)+j  k*(lx+1)*(ly+1)+(ly+1)*i+j k*(lx+1)*(ly+1)+(ly+1)*i+j+1 k*(lx+1)*(ly+1)+(ly+1)*(i-1)+j+1;];
      end
   end
end

%--------------------------------------------------------------------
% 网格显示
%--------------------------------------------------------------------

iplot=1;
if iplot==1
   figure(1)
   hold on
   axis off
   axis equal
   for ie=1:nel
        for j=1:nnel
            j1=mod(j-1,nnel)+1;
            xp(j)=x0(nodes(ie,j1),1);
            yp(j)=x0(nodes(ie,j1),2);
            zp(j)=x0(nodes(ie,j1),3);
        end
        plot3(xp,yp,zp,'-')
   end
%   pause           
end

%--------------------------------------------------------------------
%%%计算单元刚度矩阵、刚度矩阵组装
%--------------------------------------------------------------------
K=sparse(3*(ly+1)*(lx+1)*(lz+1),3*(ly+1)*(lx+1)*(lz+1));%定义刚度矩阵 
F=sparse(3*(ly+1)*(lx+1)*(lz+1),1);%定义力矩阵 
U=zeros(3*(ly+1)*(lx+1)*(lz+1),1); %定义位移矩阵 
%单元对应的节点编号
for ii0=1:nel
     
  xi1= x0(nodes(ii0,1),1);yi1=x0(nodes(ii0,1),2);zi1=x0(nodes(ii0,1),3);
  xi2=x0(nodes(ii0,2),1);yi2=x0(nodes(ii0,2),2);zi2=x0(nodes(ii0,2),3);
  xi3= x0(nodes(ii0,3),1);yi3=x0(nodes(ii0,3),2);zi3=x0(nodes(ii0,3),3);
  xi4= x0(nodes(ii0,4),1);yi4=x0(nodes(ii0,4),2);zi4=x0(nodes(ii0,4),3);
  xi5=x0(nodes(ii0,5),1);yi5=x0(nodes(ii0,5),2);zi5=x0(nodes(ii0,5),3);
  xi6= x0(nodes(ii0,6),1);yi6=x0(nodes(ii0,6),2);zi6=x0(nodes(ii0,6),3);
   xi7=x0(nodes(ii0,7),1);yi7=x0(nodes(ii0,7),2);zi7=x0(nodes(ii0,7),3);
  xi8= x0(nodes(ii0,8),1);yi8=x0(nodes(ii0,8),2);zi8=x0(nodes(ii0,8),3);
  Loc=[xi1,yi1,zi1;xi2,yi2,zi2;xi3,yi3,zi3;xi4,yi4,zi4;xi5,yi5,zi5;xi6,yi6,zi6;xi7,yi7,zi7;xi8,yi8,zi8;];
  [KE,B,D]=Stiffnesske(E,NU,Loc); %计算单元刚度
  edof=[3*nodes(ii0,1)-2,3*nodes(ii0,1)-1,3*nodes(ii0,1),3*nodes(ii0,2)-2,3*nodes(ii0,2)-1,3*nodes(ii0,2),3*nodes(ii0,3)-2,3*nodes(ii0,3)-1,3*nodes(ii0,3),3*nodes(ii0,4)-2,3*nodes(ii0,4)-1,3*nodes(ii0,4),3*nodes(ii0,5)-2,3*nodes(ii0,5)-1,3*nodes(ii0,5),3*nodes(ii0,6)-2,3*nodes(ii0,6)-1,3*nodes(ii0,6) ,3*nodes(ii0,7)-2,3*nodes(ii0,7)-1,3*nodes(ii0,7) ,3*nodes(ii0,8)-2,3*nodes(ii0,8)-1,3*nodes(ii0,8)];
 
  K(edof,edof)=K(edof,edof)+KE;%并组装刚度矩阵
end
%-------------------------------------------------------------------
%      边界处理并求解：
%------------------------------------------------------------------- 
 for i=1:lz+1
 F(3*(ly+1)*(lx+1)-1+(i-1)*3*(ly+1)*(lx+1),1) = -1; %力
 end

alldofs     = [1:3*(ly+1)*(lx+1)*(lz+1)];
for i=1:lz+1

fixeddofs   = [1+(i-1)*3*(ly+1)*(lx+1):1:3*(ly+1)+(i-1)*3*(ly+1)*(lx+1)];%固定约束
end

%-------------------------------------------------------------------
% 求解应力和应变
%-------------------------------------------------------------------
freedofs    = setdiff(alldofs,fixeddofs);
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);  %计算位移
stresplt=zeros(3*(lx+1)*(ly+1)*(lz+1),1);
for ii0=1:nel
  xi1= x0(nodes(ii0,1),1);yi1=x0(nodes(ii0,1),2);zi1=x0(nodes(ii0,1),3);
  xi2=x0(nodes(ii0,2),1);yi2=x0(nodes(ii0,2),2);zi2=x0(nodes(ii0,2),3);
  xi3= x0(nodes(ii0,3),1);yi3=x0(nodes(ii0,3),2);zi3=x0(nodes(ii0,3),3);
  xi4= x0(nodes(ii0,4),1);yi4=x0(nodes(ii0,4),2);zi4=x0(nodes(ii0,4),3);
  xi5=x0(nodes(ii0,5),1);yi5=x0(nodes(ii0,5),2);zi5=x0(nodes(ii0,5),3);
  xi6= x0(nodes(ii0,6),1);yi6=x0(nodes(ii0,6),2);zi6=x0(nodes(ii0,6),3);
   xi7=x0(nodes(ii0,7),1);yi7=x0(nodes(ii0,7),2);zi7=x0(nodes(ii0,7),3);
  xi8= x0(nodes(ii0,8),1);yi8=x0(nodes(ii0,8),2);zi8=x0(nodes(ii0,8),3);
  Loc=[xi1,yi1,zi1;xi2,yi2,zi2;xi3,yi3,zi3;xi4,yi4,zi4;xi5,yi5,zi5;xi6,yi6,zi6;xi7,yi7,zi7;xi8,yi8,zi8;];
  [KE,B,D]=Stiffnesske(E,NU,Loc); %计算应变矩阵
  edof=[3*nodes(ii0,1)-2,3*nodes(ii0,1)-1,3*nodes(ii0,1),3*nodes(ii0,2)-2,3*nodes(ii0,2)-1,3*nodes(ii0,2),3*nodes(ii0,3)-2,3*nodes(ii0,3)-1,3*nodes(ii0,3),3*nodes(ii0,4)-2,3*nodes(ii0,4)-1,3*nodes(ii0,4),3*nodes(ii0,5)-2,3*nodes(ii0,5)-1,3*nodes(ii0,5),3*nodes(ii0,6)-2,3*nodes(ii0,6)-1,3*nodes(ii0,6) ,3*nodes(ii0,7)-2,3*nodes(ii0,7)-1,3*nodes(ii0,7) ,3*nodes(ii0,8)-2,3*nodes(ii0,8)-1,3*nodes(ii0,8)];
  stresplt(edof,:)= stresplt(edof,:)+k*U(edof,:);
   stress(ii0,:)=D*B*U(edof,1);
   strain(ii0,:)=B*U(edof,1);
end
save strain.txt -ascii strain;
save stress.txt -ascii stress ;
save U.txt -ascii U ;

nnode=(ly+1)*(lx+1)*(lz+1);
fid_out=fopen('result.plt','w');
fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
fprintf(fid_out,'VARIABLES="x" "y" "z" "u" "v" "w" "sxx"  "syy" "szz"\n');
fprintf(fid_out,'ZONE T="flow-field", N= %8d,E=%8d,ET=BRICK, F=FEPOINT\n',nnode,nel);


for i=1:nnode
  fprintf(fid_out,'%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f\n', x0(i,1), x0(i,2),x0(i,3), U(3*i-2,1),  U(3*i-1,1),U(3*i,1),  stresplt(3*i-2,1),stresplt(3*i-1,1),stresplt(3*i,1));
end


for i=1:nel
    for j=1:8
     
          fprintf(fid_out,'%8d',nodes(i,j));  
    end
    fprintf(fid_out,'\n');
end;




function [KE,B,D]=Stiffnesske(E,NU,Loc)
% %
% 目的：等参单元刚度矩阵
% % 
gsx=[-0.7745966692 0 0.7745966692]; %高斯求积坐标点与系数
gsw=[0.55555555556 0.888888888889 0.55555555556];
D=[1-NU NU NU 0 0 0;NU 1-NU NU 0 0 0;NU NU 1-NU 0 0 0;0 0 0 0.5-NU 0 0;0 0 0 0 0.5-NU 0;0 0 0 0 0 0.5-NU;];
D=D*(E/((1+NU)*(1-2*NU))); %弹性矩阵
KE=zeros(24,24);
for ii=1:3 %三维高斯求积
    sx=gsx(ii);
    sw=gsw(ii);
    for jj=1:3
        nx=gsx(jj);
        nw=gsw(jj);
        for kk=1:3
            tx=gsx(kk);
            tw=gsw(kk);
            [BD,B,D]=BDcalc(sx,nx,tx,Loc,E,NU);
            KE=KE+sw*nw*tw*BD;
        end
    end
end

function [BD,B,D]=BDcalc(s,n,t,Loc,E,NU)
% % 
% B矩阵求解
% % 
dNsnt=[-(1-n)*(1-t)/8,  -(1-s)*(1-t)/8,  -(1-s)*(1-n)/8; 
    (1-n)*(1-t)/8,   -(1+s)*(1-t)/8,  -(1+s)*(1-n)/8;
    (1+n)*(1-t)/8,   (1+s)*(1-t)/8,  -(1+s)*(1+n)/8;
    -(1+n)*(1-t)/8,   (1-s)*(1-t)/8,  -(1-s)*(1+n)/8;
    -(1-n)*(1+t)/8,  -(1-s)*(1+t)/8,   (1-s)*(1-n)/8;
    (1-n)*(1+t)/8,  -(1+s)*(1+t)/8,   (1+s)*(1-n)/8;
    (1+n)*(1+t)/8,   (1+s)*(1+t)/8,   (1+s)*(1+n)/8;
    -(1+n)*(1+t)/8,  (1-s)*(1+t)/8,    (1-s)*(1+n)/8;];
dNsnt=dNsnt';
J=dNsnt*Loc;
detJ=det(J);

dNxyz=J\dNsnt;
B=zeros(6,24);
for ii=1:8 %计算B矩阵
    Bii=[dNxyz(1,ii) 0 0;0 dNxyz(2,ii) 0;0 0 dNxyz(3,ii);
        dNxyz(2,ii) dNxyz(1,ii) 0;
        0 dNxyz(3,ii) dNxyz(2,ii);
        dNxyz(3,ii) 0 dNxyz(1,ii);];
    B(:,3*(ii-1)+1:3*ii)=Bii;
end

D=[1-NU NU NU 0 0 0;NU 1-NU NU 0 0 0;NU NU 1-NU 0 0 0;0 0 0 0.5-NU 0 0;0 0 0 0 0.5-NU 0;0 0 0 0 0 0.5-NU;];
D=D*(E/((1+NU)*(1-2*NU))); %弹性矩阵

BD=detJ*transpose(B)*D*B;  %BD矩阵