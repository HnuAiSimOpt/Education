function rectangle_grid(JMP,JXYV) % 生成四面体单元
% load msh_sl
hold on;

E=length(JMP); % 读取结点信息矩阵的大小
N=length(JXYV);% 读取结点坐标矩阵的大小
JX=JXYV(:,1); % 建立总的向量
JY=JXYV(:,2);%建立总的y向量
ele_num_prnt=1; %input('是否显示单元号？1 Y or  0 N ');%是否显示单元号码
pnt_num_prnt=1; %input('是否显示结点号？1 Y  or 0 N ');%是否显示结点号码
xmin=min(JX);xmax=max(JX);x_cen=0.5*(xmin+xmax);%中点坐标
ymin=min(JY);ymax=max(JY);y_cen=0.5*(ymin+ymax);
Dx=xmax-xmin; Dy=ymax-ymin;

if Dx<Dy,xmin=x_cen-Dy/2;xmax=x_cen+Dy/2;
end
if Dx<Dy,ymin=y_cen-Dy/2;ymax=y_cen+Dy/2;
end

hold on; axis equal;

axis([xmin,xmax,ymin,ymax])
xlabel('plot of triangular grid ');hold on
del_x=0.1;del_y=0.1;%结点号放置的位置
for k=1:E
    for l=1:4%读取一个单元中的对应点的坐标
        p=JMP(k,l);
        xx(l)=JX(p);
        yy(l)=JY(p);
    end

    xx(5)=xx(1) ;
    yy(5)=yy(1);
    plot(xx,yy);
    x_cen=sum(xx(1:4))/4;
    y_cen=sum(yy(1:4))/4;
    if (ele_num_prnt==1)    %标记每个单元的单元号码
        text(x_cen-del_x,y_cen-del_y,int2str(k));
    end
end

if (pnt_num_prnt==1 )   %打印结点号
     for n=1 :N
         text(JX(n),JY(n),['(',int2str(n),')'])
     end

end
axis off