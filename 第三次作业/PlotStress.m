function PlotStress(iStress)

switch iStress

case 1

title ='x������Ӧ��';            %ͼ�α�����ʾ��x������Ӧ����

case 2

title ='y������Ӧ��';         %ͼ�α�����ʾ��y������Ӧ����

case 3

title ?='��Ӧ��';               %ͼ�α�����ʾ����Ӧ����

case 4

title ='�����Ӧ��';           %ͼ�α�����ʾ�������Ӧ����

case 5

title ='��С��Ӧ��';            %ͼ�α�����ʾ����С��Ӧ����

 end

figure;                         %����ͼ��

axis equal ;                 %����������

axis off ;                     %�ر�������

set(gcf, 'NumberTitle','off');     %�ر�NumberTitle

set(gcf,'Name',title) ;            %ͼ�α�����ʾtitle�ķ���ֵ

 for ie=1:1:16                              %����16����Ԫ����ӦӦ����ͼ

x=[gElementcoordinate(ie,1);gElementcoordinate(ie,3);gElementcoordinate(ie,5)];

y=[gElementcoordinate(ie,2);gElementcoordinate(ie,4);gElementcoordinate(ie,6)];

c=[gElementStress(ie,iStress);gElementStress(ie,iStress);gElementStress(ie,iStress)];

patch(x,y,c)

end

