%  ѧ�ţ�S230200227
%  ��������ΰ��
%  ���������ε�Ԫ��������ͨ��������
%  ��һ���߳�Ϊ2m���������м俪һ���뾶Ϊ0.4m�����Ŀס�
%  ���Ŀ��׷��壬�ϲ��ܵ����ҵĴ�СΪ200KN�ľ��ȼ��к��أ��²���֧��


% ���񻮷�
clear all;
clc;

load('Mesh.mat')


% ������Ϣ
eps=1e-10;
DirBound=@(x) (abs(x(:,2))<=eps);    % �ж�һ�����Ƿ���λ�Ʊ߽���
Dirdis=@(x) ones(size(x,1),1)*[0 0]; % λ�Ʊ߽��ϵ�λ�ƺ���
TraBound=@(x) (abs(x(:,2)-2)<=eps);  % �ж�һ�����Ƿ��ڷ��������߽���
Trafun=@(x) [1 0];                  % ��0��������

bodyforce=@(x) [0;0];          % �����
h=1;                           % ��ĺ��
E=100;nu=0.3;                  % ����ģ���Ͳ��ɱ�
D=@(x) E/(1-nu^2)*[1 nu 0;nu 1 0;0 0 (1-nu)/2];     % ƽ��Ӧ������


% �������
numn=size(node,1);
numDOF=2*numn;
Kglobal=zeros(numDOF,numDOF);
Fglobal=zeros(numDOF,1);

[Kglobal,Fglobal]=Plain_stiff(Kglobal,Fglobal,node,element,D,h,bodyforce);


% λ�Ʊ߽�����
Dirnode=find(DirBound(node));  % λ�Ʊ߽�����
DirDOF=[2*Dirnode-1;2*Dirnode];          % �߽�λ�����ɶ�
Dirdisp=reshape(Dirdis(node(Dirnode,:)),2*length(Dirnode),1);   % �߽�λ��
% ------------ �߽����  ----------------
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


% ����Ӧ���Ӧ��
[strain,stress]=postprocessing_plain(disp,D,node,element);


% ���չʾ
figure;hold on;
trimesh(element,node(:,1),node(:,2),zeros(numn,1),'edgecolor','k','linestyle','-');
scalar=1;
trimesh(element,node(:,1)+scalar*disp(1:2:end),node(:,2)+scalar*disp(2:2:end),zeros(numn,1),'edgecolor','r','facecolor','none');
daspect([1 1 1]);
axis([-1 3 0 2.5]);
box on




