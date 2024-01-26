function rectangle_grid(JMP,JXYV) % ���������嵥Ԫ
% load msh_sl
hold on;

E=length(JMP); % ��ȡ�����Ϣ����Ĵ�С
N=length(JXYV);% ��ȡ����������Ĵ�С
JX=JXYV(:,1); % �����ܵ�����
JY=JXYV(:,2);%�����ܵ�y����
ele_num_prnt=1; %input('�Ƿ���ʾ��Ԫ�ţ�1 Y or  0 N ');%�Ƿ���ʾ��Ԫ����
pnt_num_prnt=1; %input('�Ƿ���ʾ���ţ�1 Y  or 0 N ');%�Ƿ���ʾ������
xmin=min(JX);xmax=max(JX);x_cen=0.5*(xmin+xmax);%�е�����
ymin=min(JY);ymax=max(JY);y_cen=0.5*(ymin+ymax);
Dx=xmax-xmin; Dy=ymax-ymin;

if Dx<Dy,xmin=x_cen-Dy/2;xmax=x_cen+Dy/2;
end
if Dx<Dy,ymin=y_cen-Dy/2;ymax=y_cen+Dy/2;
end

hold on; axis equal;

axis([xmin,xmax,ymin,ymax])
xlabel('plot of triangular grid ');hold on
del_x=0.1;del_y=0.1;%���ŷ��õ�λ��
for k=1:E
    for l=1:4%��ȡһ����Ԫ�еĶ�Ӧ�������
        p=JMP(k,l);
        xx(l)=JX(p);
        yy(l)=JY(p);
    end

    xx(5)=xx(1) ;
    yy(5)=yy(1);
    plot(xx,yy);
    x_cen=sum(xx(1:4))/4;
    y_cen=sum(yy(1:4))/4;
    if (ele_num_prnt==1)    %���ÿ����Ԫ�ĵ�Ԫ����
        text(x_cen-del_x,y_cen-del_y,int2str(k));
    end
end

if (pnt_num_prnt==1 )   %��ӡ����
     for n=1 :N
         text(JX(n),JY(n),['(',int2str(n),')'])
     end

end
axis off