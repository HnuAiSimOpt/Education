% ���������ڱ�
% ѧ�ţ�S230200254
% ʱ�䣺20240216
close all
clear
clc

format shortEng
E = 2.000E+08;
v=0.3;
L=2;
H=0.5;
t=0.2;
scalefactor=0.2;
%�öδ��붨�������ɱ������������ϵĵ���ģ��E�����ϵĲ��ɱ�v������L������H������t����Ϊ������ֵ�����ⶨ����һ�����ڻ�ͼ�ı���scalefactor

DoE=0.125;
NEofH=round(4*H/DoE);
NEofL=round(4*L/DoE); 
nodeCoord=zeros((NEofH+1)*(NEofL+1),2);
%�滮���񻮷֣��˴��滮��������16��64=1024��

for i=1:NEofL+1
    for j=1:NEofH+1
        nodeCoord((i-1)*(NEofH+1)+j,1)=(L/NEofL)*(i-1);
        nodeCoord((i-1)*(NEofH+1)+j,2)=(H/NEofH)*(j-1);
    end
end

%nodeCoord
EleNode=zeros(NEofL*NEofH*2,3);

for i=1:NEofL
    for j=1:NEofH
        EleNode((i-1)*NEofH*2+(2*j-1),1)=(NEofH+1)*(i-1)+1+(j-1);
        EleNode((i-1)*NEofH*2+(2*j-1),2)=(NEofH+1)*i+1+(j-1);
        EleNode((i-1)*NEofH*2+(2*j-1),3)=(NEofH+1)*(i-1)+2+(j-1);
        EleNode((i-1)*NEofH*2+(2*j),1)=(NEofH+1)*(i-1)+2+(j-1);
        EleNode((i-1)*NEofH*2+(2*j),2)=(NEofH+1)*i+1+(j-1);
        EleNode((i-1)*NEofH*2+(2*j),3)=(NEofH+1)*i+2+(j-1);
    end
end
EleNode
%���δ���Խṹ��������Ԫ���񻮷֣����񻮷���ɺ����е�����Ԫ�ڵ������洢�ھ���nodeCoord�У�
%���е�Ԫ�ڵ��Ŵ洢�ھ���EleNode�У���������1105���ڵ㣬ÿ���ڵ���x��y�������꣬���nodeCoordΪ1105��2�ľ���
%������һ����16��64��2=2048�������ε�Ԫ��ÿ����Ԫ�������ڵ㣬���EleNodeΪ2048��3�ľ���
%�ڵ�͵�Ԫ���˳��Ϊ�����ң����µ���
numEle=size(EleNode,1);
numNode=size(nodeCoord,1);
%�öδ��붨������������numEle��numNode,�ֱ����ڴ洢�ṹ�еĵ�Ԫ����(�������е�EleNode����)�ͽڵ�����(������nodeCoord������)��
%���numEle=2048,numNode=1105
restrainedDof=zeros(1,2*(NEofH+1));
for i=1:2*(NEofH+1)
    restrainedDof(i)=i;
end
%restrainedDof
%�������У������ΪǶ�̣���˽�����һ�нڵ�(�ڵ�1��2��3��4��5��6��7��8��9��10��11��12��13��14��15��16��17)������ƽ�����ɶȽ���Լ����
%��Լ�������ɶȱ�Ŵ洢��restrainedDof�У�����Ϊ1��34����
xx=nodeCoord(:,1);
yy=nodeCoord(:,2);
%�öδ��붨��������������xx��yy,�ֱ����ڴ洢nodeCoord�����н���x��������нڵ��y����
numDOF=numNode*2;
%�öδ��붨����һ������numDOF�����ڴ洢���ɶȵ�����������������numNode=1105���ڵ㣬ÿ���ڵ���2�����ɶȣ����numDOF��ֵΪ1105��2=2210
displacement=zeros(numDOF,1);
%�öδ��붨����һ��������displacement,���ڴ洢����λ�ơ�����һ����numDOF=2210�����ɶȣ�����ڵ��ʼλ��Ϊ0��
%����zeros(numDof,1)������һ��numDOF��1����������������displacement
force=zeros(numDOF,1);
%�öδ��붨����һ��һ��2210��1��������force,���ڴ洢����ڵ���
P=-1000;
force(numDOF-NEofH)=P;
%force
%�öδ������������˼�����ء����Ͼ����غɡ��������غɵ�ʩ�ӣ��˴����������˼����غɵ�ʩ�ӡ�
%�����������۶�1/2���ߴ�ʩ����-z����ļ�����1000KN���ڽṹ�������ɶȾ����в����Ӧ���������ɶȱ��Ϊ(numDOF-NEofH),
%ʩ����Ӧ���غɣ�����Ԫ�ص�ֵ��Ϊ0
stiffness=zeros(numDOF);
D=(E/(1-v^2))*[1,v,0;v,1,0;0,0,(1-v)/2];
for i=1:numEle
    noindex=EleNode(i,:);
    xy=zeros(3,2);
    for j=1:3
        xy(j,1)=xx(noindex(j));
        xy(j,2)=yy(noindex(j));
    end
    J=CST_J(xy);
    Ae=1/2*det(J);
    A=1/det(J)*[J(2,2),-J(1,2),0,0;0,0,-J(2,1),J(1,1);-J(2,1),J(1,1),J(2,2),-J(1,2)];
    G=[1,0,0,0,-1,0;0,0,1,0,-1,0;0,1,0,0,0,-1;0,0,0,1,0,-1];
    B=A*G;
    eleK=t*Ae*B'*D*B;
    eleDof=[noindex(1)*2-1,noindex(1)*2,noindex(2)*2-1,noindex(2)*2,noindex(3)*2-1,noindex(3)*2];
    stiffness(eleDof,eleDof)=stiffness(eleDof,eleDof)+eleK;
end
%stiffness
%������һ����numDOF�����ɶȣ���˸öδ��붨����һ��numDOF��numDOF�ķ���stiffness,���ڴ洢�ṹ����նȾ���
%�öδ�����Ҫ����������Ҫ���裺
%1.�����е�Ԫ���б�������ø���Ԫ�ľ���ն�eleK
%2.���ݸ���Ԫ���ɶ����������ɶȵı�Ŷ�Ӧ��ϵ��������Ԫ�նȾ�����װ�ɽṹ����նȾ���
activeDof=setdiff([1:numDOF]',restrainedDof);
%�öδ��붨����һ��������activeDof�����ڴ洢���ڼ���״̬�ĵ����ɶȣ�����״̬�����ɶ�����Ϊ2210-34=2176��
%disp('Displacement')
displacement(activeDof)=stiffness(activeDof,activeDof)\force(activeDof);
%ͨ�����Ϲ��̣��Ѿ�����˽ṹ������նȾ���stiffness(eleDof,eleDof)���ڵ�������force�ͽڵ�λ������displacement,
%�öδ�����á����л��з�����ȡδԼ�����ɶȵĽڵ�λ�ƣ�ͨ���øնȾ���stiffness(activeDof,activeDof)�ҳ��ڵ�������force
%�õ��������ɶȵĽڵ�λ�ƽ������displacement(activeDof)
disp_Node=zeros(numNode,3);
for i=1:numNode
    disp_Node(i,:)=[i, displacement(i*2-1:i*2)'];
end
%disp_Node
%�öδ��밴�ڵ���������ڵ��λ�ƣ����е�һ��Ϊ�ڵ��ţ��ڶ��У�������������x��y��ƽ��λ��
stress_Node=zeros(numEle*3,5);
for i=1:numEle
    noindex=EleNode(i,:);
    xy_node=zeros(3,2);
    d=zeros(6,1);
    for j=1:3
        xy_node(j,1)=xx(noindex(j));
        xy_node(j,2)=yy(noindex(j));
        d(2*j-1)=displacement((noindex(j)-1)*2+1);
        d(2*j)=displacement((noindex(j)-1)*2+2);
    end
    J=CST_J(xy_node);
    A=1/det(J)*[J(2,2),-J(1,2),0,0;0,0,-J(2,1),J(1,1);-J(2,1),J(1,1),J(2,2),-J(1,2)];
    G=[1,0,0,0,-1,0;0,0,1,0,-1,0;0,1,0,0,0,-1;0,0,0,1,0,-1];
    B=A*G;
    sigma=D*B*d;
    stress_Node((i-1)*3+1,:)=[i,noindex(1),sigma(1),sigma(2),sigma(3)];
    stress_Node((i-1)*3+2,:)=[i,noindex(2),sigma(1),sigma(2),sigma(3)];
    stress_Node((i-1)*3+3,:)=[i,noindex(3),sigma(1),sigma(2),sigma(3)];
end
%disp('nodal stress')
%stress_Node
%�öδ������������ڵ��Ӧ�������洢��stress_Node�С����о����һ�����ڴ洢��Ԫ��ţ��ڶ������ڴ洢��Ԫ�Ľڵ㣬
%�ӵ����������зֱ����ڴ洢���ڵ��sigmax,sigmay��tauxy,�öδ�����������˸���Ԫ���ֵ�Ӧ�������
for i=1:size(restrainedDof,2)
    rownum=restrainedDof(i);
    reaction(i)=stiffness(rownum,:)*displacement-force(rownum);
end
%disp('Reaction')
%reaction'
%�öδ���������ȡ֧������������������6��֧��������ͨ���������ɶȶ�Ӧ������նȾ���stiffness(rownum,:)������λ������displacement��˵õ���
maxlen=-1;
for i=1:numEle
    noindex=EleNode(i,:);
    deltax=xx(noindex(2))-xx(noindex(1));
    deltay=yy(noindex(2))-yy(noindex(1));
    L=sqrt(deltax*deltax+deltay*deltay);
    if(L>maxlen)
        maxlen=L;
    end
end
maxabsDisp=0;
for i=1:numNode
    tempdispU=displacement(i*2-1);
    tempdispV=displacement(i*2);
    tempdisp=sqrt(tempdispU*tempdispU+tempdispV*tempdispV);
    if(tempdisp>maxabsDisp)
        maxabsDisp=tempdisp;
    end
end
factor=0.1;
if(maxabsDisp>1e-30)
    factor=scalefactor*maxlen/maxabsDisp;
end
figure
for i=1:numEle
    noindex1=EleNode(i,:);
    noindex=[noindex1,noindex1(1)];
    Line1=plot(xx([noindex]),yy([noindex]),'--','Color',[0.4,0.4,0.4],'LineWidth',2);
    hold on
end
xxDeformed=xx;
yyDeformed=yy;
for i=1:numNode
    xxDeformed(i)=xxDeformed(i)+factor*displacement(i*2-1);
    yyDeformed(i)=yyDeformed(i)+factor*displacement(i*2);
end
for i=1:numNode
    xxDeformed(i)=xxDeformed(i)+factor*displacement(i*2-1);
    yyDeformed(i)=yyDeformed(i)+factor*displacement(i*2);
end
for i=1:numEle
    noindex=EleNode(i,:);
    coordx=[xxDeformed(noindex(1)),xxDeformed(noindex(2)),xxDeformed(noindex(3)),xxDeformed(noindex(1))];
    coordy=[yyDeformed(noindex(1)),yyDeformed(noindex(2)),yyDeformed(noindex(3)),yyDeformed(noindex(1))];
    stressX=[stress_Node((i-1)*3+1,3),stress_Node((i-1)*3+2,3),stress_Node((i-1)*3+3,3),stress_Node((i-1)*3+3,4)];
    fill(coordx,coordy,stressX);
     shading interp;
    colorbar;
end
 for i=1:numEle
     noindex1=EleNode(i,:);
     noindex=[noindex1,noindex1(1)];
     Line2=plot(xxDeformed([noindex]),yyDeformed([noindex]),'-','LineWidth',0.5);
     hold on
 end
minx=-0.2*(max(xxDeformed)-min(xxDeformed))+min(xxDeformed);
maxx=0.2*(max(xxDeformed)-min(xxDeformed))+max(xxDeformed);
miny=-0.2*(max(yyDeformed)-min(yyDeformed))+min(yyDeformed);
maxy=0.2*(max(yyDeformed)-min(yyDeformed))+max(yyDeformed);
axis([minx,maxx,miny,maxy]);
xlabel('X');
ylabel('Y');
zlabel('Z');
legend([Line1,Line2],'Undeformed Shape','Deformed Shape')
%�öδ������ڻ��ƽṹ�ı���ͼ����Ҫ������������:
%1.��ͼ��������:�öδ����������ø����ͼ��������Ҫ����:��ø�����Ԫ����󳤶ȵ�maxlen;��ýڵ�λ�����ֵmaxabsDisp;
%����������������õĵ�Ԫ��󳤶Ⱥͽڵ�λ�Ƶ����ֵ����������ͼ�λ��Ƶ�����ϵ��factor
%2.���Ʊ���ͼ:�öδ������ڻ��ƽṹ����ͼ����Ҫ������������:�������е�Ԫ�����ݵ�Ԫ��ڵ�Ĺ�ϵ��������ߣ����Ʊ���֮ǰ�Ľṹ;
%���ݽṹ���κ�ĵĽڵ����꣬���Ʊ��κ�Ľṹ��״;�ڸ��ڵ㴦������Ӧ�Ľڵ���;�ڸ���Ԫ��������Ӧ�ĵ�Ԫ���;ͨ������x��y����ʾ��Χ��
%ʹ��ͼ���Ժ��ʵı���������ʾ
%3.���ƽڵ�Ӧ��:�öδ��������ڱ���ͼ�Ļ����ϻ��ƽڵ�Ӧ������Ҫ������������:���ƽڵ�Ӧ��sigmax;���ƽڵ�Ӧ��sigmay;���ƽڵ�Ӧ��tauxy