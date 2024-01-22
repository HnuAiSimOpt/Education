% 问题描述
% 平面弹性问题，一块板的一条边固定，向另一个顶点施加平行于固定边的力。
% 采用线性三角单元。
% by S230200235 雷沐昆

%% 网格划分
a=2;b=1;p=-1;
Ny=5;Nx=2*Ny;
[node,element]=Rect_mesh(a,b,Nx,Ny,'tri');

%% 问题信息
DirDOF=(1:2*(Ny+1))';          % 位移自由度
Dirdisp=zeros(size(DirDOF));   % 边界位移

bodyforce=@(x) [0;0];          % 体积力
h=1;                           % 板的厚度
E=4;nu=0;
D=@(x) E*(1-nu)/(1+nu)/(1-2*nu)*[1 nu/(1-nu) 0;nu/(1-nu) 1 0;0 0 (1-2*nu)/2/(1-nu)];                     % 平面应变问题
% D=@(x) E/(1-nu^2)*[1 nu 0;nu 1 0;0 0 (1-nu)/2];                     % 平面应力问题

NeuDOF=2*((Ny+1)*Nx+(1:Ny+1));         % 非零的力边界自由度
NeuF=b/Ny*p*h*ones(length(NeuDOF),1);  % 结点力
NeuF(1)=NeuF(1)/2;NeuF(end)=NeuF(end)/2;

%% 整体分析
numnode=size(node,1);
numDOF=2*numnode;
Kglobal=zeros(numDOF,numDOF);
Fglobal=zeros(numDOF,1);

[Kglobal,Fglobal]=plane_assembly_tri(Kglobal,Fglobal,node,element,bodyforce,D,h);


% ------------ 边界荷载  ----------------
Fglobal(NeuDOF)=Fglobal(NeuDOF)+NeuF;
% ------------------------------------------
%% 考虑位移边界条件
AllDOF=1:numDOF;
EffDOF=truncat(AllDOF,DirDOF);
Keff =Kglobal(EffDOF,EffDOF);
Kd=Kglobal(EffDOF,DirDOF);
Feff=Fglobal(EffDOF)-Kd*Dirdisp;
deff=Keff\Feff;
disp=zeros(numDOF,1);
disp(DirDOF)=Dirdisp;
disp(EffDOF)=deff;

%% 计算应变和应力
[strain,stress]=postprocessing_plain(disp,D,node,element);

strain_C=disp(end-1:end)';
stress_C=stress(end,:);
%% 结果展示
figure;hold on;
trimesh(element,node(:,1),node(:,2),zeros(numnode,1),'edgecolor','k','linestyle','--');
scalar=0.01;
trimesh(element,node(:,1)+scalar*disp(1:2:end),node(:,2)+scalar*disp(2:2:end),zeros(numnode,1),'edgecolor','r','facecolor','none');

% end



