%������Ԭ�ұ�  ѧ�ţ�S230200228 
%purpose�������������ε�Ԫ�Ա���ṹ����Ӧ��Ӧ�����
%�߽�������λ�Ʊ߽��������ڵ�1��4λ��Ϊ0   ���߽��������ڵ�2��3��y��������Ϊ0��
clc
clear;
numele=2;%��Ԫ��
numnode=4;%��Ԫ�ڵ���
elenodecorre=[1 3 4;1 2 3];%������Ԫ�ֱ��Ӧ�Ľڵ�
%�Ƚ����غ�ǰ��ͼ������
ininodevector=[0 0 0.5 0 0.5 0.25 0 0.25];%�����ڵ������
figure
axis off
axis equal
hold on
for inicoori=1:numnode
    plot(ininodevector(2*inicoori-1),ininodevector(2*inicoori),'bo');
end
for inicoorj=1:numele
    hold on
    line([ininodevector(2*elenodecorre(inicoorj,1)-1) ...
        ininodevector(2*elenodecorre(inicoorj,2)-1)],...
        [ininodevector(elenodecorre(inicoorj,1)*2) ...
        ininodevector(elenodecorre(inicoorj,2)*2)]);
    hold on
    line([ininodevector(2*elenodecorre(inicoorj,2)-1) ...
        ininodevector(2*elenodecorre(inicoorj,3)-1)],...
        [ininodevector(elenodecorre(inicoorj,2)*2) ...
        ininodevector(elenodecorre(inicoorj,3)*2)]);
    hold on
    line([ininodevector(2*elenodecorre(inicoorj,3)-1) ...
        ininodevector(2*elenodecorre(inicoorj,1)-1)],...
        [ininodevector(elenodecorre(inicoorj,3)*2) ...
        ininodevector(elenodecorre(inicoorj,1)*2)]);
end
%�ȸ�����ز���
E=210e6;%����ģ��
NU=0.3;%���ɱ�
t=0.025;%�����ȣ�
k1=LinearTriangleElementStiffness(E,NU,t,0,0,0.5,0.25,0,0.25,1);%�ֱ�������������ε�Ԫ�ĸնȾ���k1��k2
k2=LinearTriangleElementStiffness(E,NU,t,0,0,0.5,0,0.5,0.25,1);
K=zeros(8,8);%��Ϊ���ĸ��ڵ㣬����ĸնȾ���Ϊ8x8������������һ��8x8�Ŀվ���
K=LinearTriangleAssemble(K,k1,1,3,4);
K=LinearTriangleAssemble(K,k2,1,2,3);
%������Ҫ���б߽���������
%U1x=U1y=U4x=U4y=0----λ�Ʊ߽�������1��4�ڵ�λ��Ϊ0
%F2x=F3x=9.375;F2y=F3y=0----���߽�������
k=K(3:6,3:6);
f=[9.375;0;9.375;0];
u=k\f;%�ø�˹��ȥ�������2��3λ��
%�õ��ڵ�2��3λ�ƺ����ڵ�1��4λ�Ʒ���
U=[0;0;u;0;0];
F=K*U;
u1=[U(1);U(2);U(5);U(6);U(7);U(8)];
u2=[U(1);U(2);U(3);U(4);U(5);U(6)];
sigma1=LinearTriangleElementStresses(E,NU,t,0,0,0.5,0.25,0,0.25,1,u1);
sigma2=LinearTriangleElementStresses(E,NU,t,0,0,0.5,0,0.5,0.25,1,u2);%�õ���Ԫ1��2�Ľڵ�Ӧ��ֵ
s1=LinearTriangleElementPStresses(sigma1);
s2=LinearTriangleElementPStresses(sigma2);
%��صļ����������ұߵĶ�Ӧ�����
finnodevector=[0 0 0.5+u(1,1) 0+u(2,1) 0.5+u(3,1) 0.25-u(4,1) 0 0.25];%�ܵ��غɸ����ڵ������
hold on
%�������غɺ��ͼ�����ڱ仯̫С��ͼ����ʾ������
for fincoori=1:numnode
    plot(finnodevector(2*fincoori-1),finnodevector(2*fincoori),'ro');
end
for fincoorj=1:numele
    hold on
    line([finnodevector(2*elenodecorre(fincoorj,1)-1) ...
        finnodevector(2*elenodecorre(fincoorj,2)-1)],...
        [finnodevector(elenodecorre(fincoorj,1)*2) ...
        finnodevector(elenodecorre(fincoorj,2)*2)],'Color','red','LineStyle','--');
    hold on
    line([finnodevector(2*elenodecorre(fincoorj,2)-1) ...
        finnodevector(2*elenodecorre(fincoorj,3)-1)],...
        [finnodevector(elenodecorre(fincoorj,2)*2) ...
        finnodevector(elenodecorre(fincoorj,3)*2)],'Color','red','LineStyle','--');
    hold on
    line([finnodevector(2*elenodecorre(fincoorj,3)-1) ...
        finnodevector(2*elenodecorre(fincoorj,1)-1)],...
        [finnodevector(elenodecorre(fincoorj,3)*2) ...
        finnodevector(elenodecorre(fincoorj,1)*2)],'Color','red','LineStyle','--');
end