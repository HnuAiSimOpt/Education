function plotdisp(nel,disp,nodes,x0,idisp) % ��ʾӦ����ͼ
% Nz=1037
Nz=1037

% iStress Ӧ������ָ�꣬�������������ֵ
% 1����u����λ��
% 2����v����λ��  

%% ȷ��ͼ�δ�������
switch idisp
    case 1
        title='x�����ٶ���ͼ';
    case 2
        title='v����λ��';
end
u=zeros(Nz,2);
for i=1:Nz
    u(i,1)=disp(i);
%     u(i,2)=disp(2*i);
end

%% ����ͼ�δ��ڣ������������ᣬָ������
figure;
axis equal;
axis off;
set(gcf,'NumberTitle','off');
set(gcf,'Name',title);               
uMin=min(u(:,idisp));  %ȷ����Сλ��
uMax=max(u(:,idisp));  %���λ��
caxis([uMin,uMax]);           %ָ����ɫӳ���
colormap('jet');

%% ���ݵ�Ԫ�ڵ������Ӧ��ֵ���������ı��Σ���ʾӦ����ͼ
for ie=1:1:nel
    x=[x0(nodes(ie,1),1);
       x0(nodes(ie,2),1);
       x0(nodes(ie,3),1);
       x0(nodes(ie,4),1)];
    y=[x0(nodes(ie,1),2);
       x0(nodes(ie,2),2);
       x0(nodes(ie,3),2);
       x0(nodes(ie,4),2)];
   c=[ u(nodes(ie,1),idisp);
       u(nodes(ie,2),idisp);
       u(nodes(ie,3),idisp);
       u(nodes(ie,4),idisp)];
   set(patch(x,y,c),'EdgeColor','interp')
end

%% ������ɫ��������ָʾ��ͬ��ɫ����Ӧ��Ӧ��ֵ
yTick=uMin:(uMax-uMin)/10:uMax;
Label=cell(1,length(yTick));
for i=1:length(yTick)
    Label{i}=sprintf('%2f',yTick(i));
end
set(colorbar('vert'),'YTick',yTick,'YTickLabelMode','Manual','YTickLabel',Label);