% 六面体单元
%----------------------------------------------------------------------
%%  主要参数
% goords  节点坐标值
% nodes   单元节点编号
% nel     单元总量
% KE      单元刚度矩阵
% K       系统刚度矩阵
% F       系统载荷向量
% disp    系统节点位移向量
% stress  系统应力矩阵
% strain  应变矩阵
%-----------------------------------------------------------------------


%%  前处理

%获取模型节点坐标值及坐标编号
[goords,nodes]=getcoord();
 % 总单元数
nel=length(nodes(:,1));
%总节点数
nnode=length(goords(:,1))
%系统总自由度
sdof=3*length(goords(:,1));
%单元自由度
edof=24;
%初始化刚度、位移、载荷矩阵
KE=[edof,edof];       %单元刚度矩阵
K=sparse(sdof,sdof);  %定义系统刚度矩阵 
F=sparse(sdof,1);     %定义力矩阵 
disp=zeros(sdof,1);   %定义位移矩阵 l=8; 每个元素的节点数   

%%  手动设立条件

%定义材料
E=input('请输入杨氏模量E=\n');
P=input('请输入泊松比P=\n');

%定义集中载荷数目
i=input('请输入需要施加的载荷数目=\n');

%循环读取载荷位置及大小，为该点自由度编号
for  j=1:i
    
  b=input('请输入需要施加的载荷编号=\n');
  F(b,1)=input('请输入需要施加的载荷大小=\n');
  
end


%%   刚度矩阵的建立

for j=1:nel  
   % 一次取一个单元的编号
   dex=[];
   dex=[nodes(j,2:9)] ;
   
  % 该单元节点坐标
   XYZ=[];
   for k=1:8      
     XYZ(k,1)=goords(dex(1,k),2);
      XYZ(k,2)=goords(dex(1,k),3);
       XYZ(k,3)=goords(dex(1,k),4);      
   end
   [KE,B,D]=Stiffnesske(E,P,XYZ); %计算应变矩阵
   
   index=[];
     for i=1:8      
       index=[index 3*dex(1,i)-2 3*dex(1,i)-1 3*dex(1,i)];     
     end
     
      K(index,index)=K(index,index)+KE;%并组装刚度矩阵
end

%%  边界条件及求解模块
%输入边界条件
alldofs   = [1:sdof];
bcval=[169,170,171,172,173,174,175,176,177,178,179,180,181,182 ];
c=length(bcval);
fixeddofs=[];
for i=1:c    
 fixeddofs=[fixeddofs 3*bcval(i)-2,3*bcval(i)-1,3*bcval(i)]; 
end
 n=length(fixeddofs);
% 根据置一法修改刚度矩阵和载荷矩阵
 for i=1:n
    c= fixeddofs(i);
    for j=1:sdof
       K(c,j)=0;
    end
    K(c,c)=1;
    F(c)=0;
 end
 
 %求解位移矩阵
disp= K \ F;  



%% 后处理模块
%创建“456.plt”文件并将输出结果寄存再该文件中
fid_out=fopen('result .plt','w');
fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
fprintf(fid_out,'VARIABLES="x" "y" "z" "u" "v" "w" \n');
fprintf(fid_out,'ZONE T="flow-field", N= %8d,E=%8d,ET=BRICK, F=FEPOINT\n',nnode,nel);
for i=1:nnode
  fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n', goords(i,2), goords(i,3),goords(i,4), disp(3*i-2,1)+0, disp(3*i-1,1)+0,disp(3*i,1)+0);
end
for i=1:nel
    fprintf(fid_out,'%8d%8d%8d%8d%8d%8d%8d%8d\n',nodes(i,2),nodes(i,3),nodes(i,4),nodes(i,5),nodes(i,6),nodes(i,7),nodes(i,8),nodes(i,9));
end



