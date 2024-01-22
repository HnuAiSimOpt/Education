%  学号：S230200227
%  姓名：杜伟栋
%  采用三角形单元分析的普通弹性问题
%  在一个边长为2m的正方形中间开一个半径为0.4m的中心孔。
%  中心开孔方板，上部受到向右的大小为200KN的均匀剪切荷载，下部固支。


% 网格划分
clear all;
clc;

load('Mesh.mat')


% 问题信息
eps=1e-10;
DirBound=@(x) (abs(x(:,2))<=eps);    % 判断一个点是否在位移边界上
Dirdis=@(x) ones(size(x,1),1)*[0 0]; % 位移边界上的位移函数
TraBound=@(x) (abs(x(:,2)-2)<=eps);  % 判断一个点是否在非零面力边界上
Trafun=@(x) [1 0];                  % 非0面力函数

bodyforce=@(x) [0;0];          % 体积力
h=1;                           % 板的厚度
E=100;nu=0.3;                  % 弹性模量和泊松比
D=@(x) E/(1-nu^2)*[1 nu 0;nu 1 0;0 0 (1-nu)/2];     % 平面应力问题


% 整体分析
numn=size(node,1);
numDOF=2*numn;
Kglobal=zeros(numDOF,numDOF);
Fglobal=zeros(numDOF,1);

[Kglobal,Fglobal]=Plain_stiff(Kglobal,Fglobal,node,element,D,h,bodyforce);


% 位移边界条件
Dirnode=find(DirBound(node));  % 位移边界结点编号
DirDOF=[2*Dirnode-1;2*Dirnode];          % 边界位移自由度
Dirdisp=reshape(Dirdis(node(Dirnode,:)),2*length(Dirnode),1);   % 边界位移
% ------------ 边界荷载  ----------------
[NeuDOF,NeuF]=Traction_processing(node,element,TraBound,Trafun,h);
Fglobal(NeuDOF)=Fglobal(NeuDOF)+NeuF;
% ------------------------------------------

AllDOF=1:numDOF;
EffDOF=truncat(AllDOF,DirDOF);
Keff =Kglobal(EffDOF,EffDOF);
Kd=Kglobal(EffDOF,DirDOF);
Feff=Fglobal(EffDOF)-Kd*Dirdisp;
deff=Keff\Feff;
disp=zeros(numDOF,1);
disp(DirDOF)=Dirdisp;
disp(EffDOF)=deff;


% 计算应变和应力
[strain,stress]=postprocessing_plain(disp,D,node,element);


% 结果展示
figure;hold on;
trimesh(element,node(:,1),node(:,2),zeros(numn,1),'edgecolor','k','linestyle','-');
scalar=1;
trimesh(element,node(:,1)+scalar*disp(1:2:end),node(:,2)+scalar*disp(2:2:end),zeros(numn,1),'edgecolor','r','facecolor','none');
daspect([1 1 1]);
axis([-1 3 0 2.5]);
box on




