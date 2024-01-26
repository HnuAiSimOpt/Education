clear
E = 1e7;
NU = 1/3;
t = 0.1;
ID = 1;
node = [2 1
		  2 0
		  0 1
		  0 0];
element = [3 2 1
			  2 3 4];
force = [1 2 -5000
			2 2 -5000];
constrain = [3 1 0
				 3 2 0
				 4 1 0
				 4 2 0];

%ȷ���ڵ㡢��Ԫ����
[nnode,ntmp]=size(node);
[nelem,etmp]=size(element);
[nforce,ftmp]=size(force);
[nconstrain,ctmp]=size(constrain);


KKG=zeros(2*nnode);
FFG=zeros(2*nnode,1);
UUG=zeros(2*nnode,1);
StressElem=zeros(nelem,3);
StressNode=zeros(nnode,3);
k=zeros(6,6);     % ��Ԫ�նȾ���

%��Ӧ���ϼ��������
E=1e6;
NU=1/3;
t=1;
ID=1;% ƽ��״̬��־

for i = 1:nelem
	% node(element(i,:),1) ��ʾ��i����Ԫ�Ľڵ�Ŷ�Ӧ������x����ֵ
	k=T_Stiffness(E,NU,t,node(element(i,:),1),node(element(i,:),2),ID);
	% element(i,1)��ʾ��i����Ԫ�ĵ�һ������
	KKG = Triangle2D3Node_Assembly(KKG,k,element(i,1),element(i,2),element(i,3));
end

% �߽������Ĵ����նȷ������
kk=KKG(1:4,1:4);
p=[0;-5000;0;-5000];
u=kk\p
% �ڵ���RF�ļ���
U=[u;0;0;0;0];
RF=KKG*U
% ����Ԫ��Ӧ������
for i=1:nelem
    l=element(i,1);m=element(i,2);n=element(i,3);
    ele_u=[U(2*l-1),U(2*l),U(2*m-1),U(2*m),U(2*n-1),U(2*n)]';
    stress=T_Stress(E,NU,node(element(i,:),1),node(element(i,:),2),ele_u,ID);
    StressElem(i,1)=stress(1,1);
    StressElem(i,2)=stress(2,1);
    StressElem(i,3)=stress(3,1);
end
