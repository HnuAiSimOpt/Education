% ��������
% ƽ�浯�����⣬һ����һ���߹̶�������һ������ʩ��ƽ���ڹ̶��ߵ�����
% �����������ǵ�Ԫ��
% by S230200235 ������

%% ���񻮷�
a=2;b=1;p=-1;
Ny=5;Nx=2*Ny;
[node,element]=Rect_mesh(a,b,Nx,Ny,'tri');

%% ������Ϣ
DirDOF=(1:2*(Ny+1))';          % λ�����ɶ�
Dirdisp=zeros(size(DirDOF));   % �߽�λ��

bodyforce=@(x) [0;0];          % �����
h=1;                           % ��ĺ��
E=4;nu=0;
D=@(x) E*(1-nu)/(1+nu)/(1-2*nu)*[1 nu/(1-nu) 0;nu/(1-nu) 1 0;0 0 (1-2*nu)/2/(1-nu)];                     % ƽ��Ӧ������
% D=@(x) E/(1-nu^2)*[1 nu 0;nu 1 0;0 0 (1-nu)/2];                     % ƽ��Ӧ������

NeuDOF=2*((Ny+1)*Nx+(1:Ny+1));         % ��������߽����ɶ�
NeuF=b/Ny*p*h*ones(length(NeuDOF),1);  % �����
NeuF(1)=NeuF(1)/2;NeuF(end)=NeuF(end)/2;

%% �������
numnode=size(node,1);
numDOF=2*numnode;
Kglobal=zeros(numDOF,numDOF);
Fglobal=zeros(numDOF,1);

[Kglobal,Fglobal]=plane_assembly_tri(Kglobal,Fglobal,node,element,bodyforce,D,h);


% ------------ �߽����  ----------------
Fglobal(NeuDOF)=Fglobal(NeuDOF)+NeuF;
% ------------------------------------------
%% ����λ�Ʊ߽�����
AllDOF=1:numDOF;
EffDOF=truncat(AllDOF,DirDOF);
Keff =Kglobal(EffDOF,EffDOF);
Kd=Kglobal(EffDOF,DirDOF);
Feff=Fglobal(EffDOF)-Kd*Dirdisp;
deff=Keff\Feff;
disp=zeros(numDOF,1);
disp(DirDOF)=Dirdisp;
disp(EffDOF)=deff;

%% ����Ӧ���Ӧ��
[strain,stress]=postprocessing_plain(disp,D,node,element);

strain_C=disp(end-1:end)';
stress_C=stress(end,:);
%% ���չʾ
figure;hold on;
trimesh(element,node(:,1),node(:,2),zeros(numnode,1),'edgecolor','k','linestyle','--');
scalar=0.01;
trimesh(element,node(:,1)+scalar*disp(1:2:end),node(:,2)+scalar*disp(2:2:end),zeros(numnode,1),'edgecolor','r','facecolor','none');

% end



