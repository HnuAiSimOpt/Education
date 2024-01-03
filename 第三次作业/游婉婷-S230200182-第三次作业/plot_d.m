function []=plot_d(d,node)
d1=[];
d2=[];
for i=1:222
 d1=[d1;d(2*i-1)];
 d2=[d2;d(2*i)];
end
d11=d1;  %得到水平方向各节点位移
d22=d2;  %得到竖直方向各节点位移
X1=linspace(min(node(:,2)),max(node(:,2)),200);  %将横纵坐标均分成 200 份
Y1=linspace(min(node(:,3)),max(node(:,3)),200);
[X,Y,d1]=griddata(node(:,2),node(:,3),d1,X1',Y1,'natural');  %对二维或三维散点数据插值,插值方法为'natural'
[X,Y,d2]=griddata(node(:,2),node(:,3),d2,X1',Y1,'v4');       %插值方法为'v4'
fig_1 = figure('Name','figure_1','Position',[50,500,500,500]);
subplot(4,1,1);
contourf(X,Y,d1);
colorbar  %分布柱图
title('水平方向位移分布等高图/m')
box on
subplot(4,1,2);
pcolor(X,Y,d1);
shading interp  %使色彩平滑过渡
colorbar  %分布柱图 
axis equal  %将横轴纵轴的定标系数设成相同值
title('水平方向位移分布云图/m')
box on;
x1_row=find(abs(node(:,2)-2)<1e-5);  %找到右上角和右下角的节点标号
y1_row=find(abs(node(:,3)-1)<1e-5);
y2_row=find(abs(node(:,3))<1e-5);
node_num1=intersect(x1_row,y1_row);
node_num2=intersect(x1_row,y2_row);
x2_row=find(abs(node(:,2)-0)<1e-5);  %找到左上角和左下角的节点标号
y3_row=find(abs(node(:,3)-1)<1e-5);
y4_row=find(abs(node(:,3))<1e-5);
node_num3=intersect(x2_row,y3_row);
node_num4=intersect(x2_row,y4_row);
text_x2y1=d11(node_num1);
text_x2y2=d11(node_num2);
text_x2y3=d22(node_num3);
text_x2y4=d22(node_num4);
text_1=['坐标','(',num2str(node(node_num1,2)),',',num2str(node(node_num1,3)),')',' 位移为',num2str(text_x2y1)];
text_2=['坐标','(',num2str(node(node_num2,2)),',',num2str(node(node_num2,3)),')',' 位移为',num2str(text_x2y2)];
text_3=['坐标','(',num2str(node(node_num3,2)),',',num2str(node(node_num3,3)),')',' 位移为',num2str(text_x2y3)];
text_4=['坐标','(',num2str(node(node_num4,2)),',',num2str(node(node_num4,3)),')',' 位移为',num2str(text_x2y4)];
text(node(node_num1,2),node(node_num1,3),text_1,'Color','magenta','FontSize',10);  %在云图中标注周边四个节点
text(node(node_num2,2),node(node_num2,3),text_2,'Color','magenta','FontSize',10);
text(node(node_num3,2),node(node_num3,3),text_3,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
text(node(node_num4,2),node(node_num4,3),text_4,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
subplot(4,1,3);
contourf(X,Y,d2);
colorbar  %分布柱图
title('竖直方向位移分布等高图/m')
box on;
subplot(4,1,4);
pcolor(X,Y,d2);
shading interp  %使色彩平滑过渡
colorbar  %分布柱图 
axis equal  %将横轴纵轴的定标系数设成相同值
title('竖直方向位移分布云图/m')
box on;
x1_row=find(abs(node(:,2)-2)<1e-5);  %找到右上角和右下角的节点标号
y1_row=find(abs(node(:,3)-1)<1e-5);
y2_row=find(abs(node(:,3))<1e-5);
node_num1=intersect(x1_row,y1_row);
node_num2=intersect(x1_row,y2_row);
x2_row=find(abs(node(:,2)-0)<1e-5);  %找到左上角和左下角的节点标号
y3_row=find(abs(node(:,3)-1)<1e-5);
y4_row=find(abs(node(:,3))<1e-5);
node_num3=intersect(x2_row,y3_row);
node_num4=intersect(x2_row,y4_row);
text_x2y1=d22(node_num1);
text_x2y2=d22(node_num2);
text_x2y3=d22(node_num3);
text_x2y4=d22(node_num4);
text_1=['坐标','(',num2str(node(node_num1,2)),',',num2str(node(node_num1,3)),')',' 位移为',num2str(text_x2y1)]; 
text_2=['坐标','(',num2str(node(node_num2,2)),',',num2str(node(node_num2,3)),')',' 位移为',num2str(text_x2y2)];
text_3=['坐标','(',num2str(node(node_num3,2)),',',num2str(node(node_num3,3)),')',' 位移为',num2str(text_x2y3)];
text_4=['坐标','(',num2str(node(node_num4,2)),',',num2str(node(node_num4,3)),')',' 位移为',num2str(text_x2y4)];
text(node(node_num1,2),node(node_num1,3),text_1,'Color','magenta','FontSize',10);  %在云图中标注周边四个节点
text(node(node_num2,2),node(node_num2,3),text_2,'Color','magenta','FontSize',10);
text(node(node_num3,2),node(node_num3,3),text_3,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
text(node(node_num4,2),node(node_num4,3),text_4,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
sgtitle('虚拟空心结构位移分布图/m') ;
fig_2 = figure('Name','figure_2','Position',[600,500,500,500]);
X=X';
for i=1:200
    for j=1:200
    if (X(i,j)>0.8&&X(i,j)<1.2)&&(Y(i,j)>0.3&&Y(i,j)<0.7)
        X(i,j)=NaN;
        Y(i,j)=NaN;
    end
    end
end
for i=81:120
    for j=1:80
        X(i,j)=X(i-1,j)+0.0101;
    end
end
for i=81:120
    for j=121:200
        X(i,j)=X(i-1,j)+0.0101;
    end
end
X=X';
for i=81:120
    for j=1:80
        Y(i,j)=Y(i-1,j)+0.00505;
    end
end
for i=81:120
    for j=121:200
        Y(i,j)=Y(i-1,j)+0.00505;
    end
end
subplot(4,1,1);
contourf(X,Y,d1);
colorbar  %分布柱图
title('水平方向位移分布等高图/m')
box on
subplot(4,1,2);
pcolor(X,Y,d1);
shading interp  %使色彩平滑过渡
colorbar  %分布柱图 
axis equal  %将横轴纵轴的定标系数设成相同值
title('水平方向位移分布云图/m')
box on;
x1_row=find(abs(node(:,2)-2)<1e-5);  %找到右上角和右下角的节点标号
y1_row=find(abs(node(:,3)-1)<1e-5);
y2_row=find(abs(node(:,3))<1e-5);
node_num1=intersect(x1_row,y1_row);
node_num2=intersect(x1_row,y2_row);
x2_row=find(abs(node(:,2)-0)<1e-5);  %找到左上角和左下角的节点标号
y3_row=find(abs(node(:,3)-1)<1e-5);
y4_row=find(abs(node(:,3))<1e-5);
node_num3=intersect(x2_row,y3_row);
node_num4=intersect(x2_row,y4_row);
text_x2y1=d11(node_num1);
text_x2y2=d11(node_num2);
text_x2y3=d22(node_num3);
text_x2y4=d22(node_num4);
text_1=['坐标','(',num2str(node(node_num1,2)),',',num2str(node(node_num1,3)),')',' 位移为',num2str(text_x2y1)];
text_2=['坐标','(',num2str(node(node_num2,2)),',',num2str(node(node_num2,3)),')',' 位移为',num2str(text_x2y2)];
text_3=['坐标','(',num2str(node(node_num3,2)),',',num2str(node(node_num3,3)),')',' 位移为',num2str(text_x2y3)];
text_4=['坐标','(',num2str(node(node_num4,2)),',',num2str(node(node_num4,3)),')',' 位移为',num2str(text_x2y4)];
text(node(node_num1,2),node(node_num1,3),text_1,'Color','magenta','FontSize',10);  %在云图中标注周边四个节点
text(node(node_num2,2),node(node_num2,3),text_2,'Color','magenta','FontSize',10);
text(node(node_num3,2),node(node_num3,3),text_3,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
text(node(node_num4,2),node(node_num4,3),text_4,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
subplot(4,1,3);
contourf(X,Y,d2);
colorbar  %分布柱图
title('竖直方向位移分布等高图/m')
box on;
subplot(4,1,4);
pcolor(X,Y,d2);
shading interp  %使色彩平滑过渡
colorbar  %分布柱图 
axis equal  %将横轴纵轴的定标系数设成相同值
title('竖直方向位移分布云图/m')
box on;
x1_row=find(abs(node(:,2)-2)<1e-5);  %找到右上角和右下角的节点标号
y1_row=find(abs(node(:,3)-1)<1e-5);
y2_row=find(abs(node(:,3))<1e-5);
node_num1=intersect(x1_row,y1_row);
node_num2=intersect(x1_row,y2_row);
x2_row=find(abs(node(:,2)-0)<1e-5);  %找到左上角和左下角的节点标号
y3_row=find(abs(node(:,3)-1)<1e-5);
y4_row=find(abs(node(:,3))<1e-5);
node_num3=intersect(x2_row,y3_row);
node_num4=intersect(x2_row,y4_row);
text_x2y1=d22(node_num1);
text_x2y2=d22(node_num2);
text_x2y3=d22(node_num3);
text_x2y4=d22(node_num4);
text_1=['坐标','(',num2str(node(node_num1,2)),',',num2str(node(node_num1,3)),')',' 位移为',num2str(text_x2y1)]; 
text_2=['坐标','(',num2str(node(node_num2,2)),',',num2str(node(node_num2,3)),')',' 位移为',num2str(text_x2y2)];
text_3=['坐标','(',num2str(node(node_num3,2)),',',num2str(node(node_num3,3)),')',' 位移为',num2str(text_x2y3)];
text_4=['坐标','(',num2str(node(node_num4,2)),',',num2str(node(node_num4,3)),')',' 位移为',num2str(text_x2y4)];
text(node(node_num1,2),node(node_num1,3),text_1,'Color','magenta','FontSize',10);  %在云图中标注周边四个节点
text(node(node_num2,2),node(node_num2,3),text_2,'Color','magenta','FontSize',10);
text(node(node_num3,2),node(node_num3,3),text_3,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
text(node(node_num4,2),node(node_num4,3),text_4,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
sgtitle('实体空心结构位移分布图/m');
fig_3 = figure('Name','figure_3','Position',[1150,500,500,500]);
node_num=0;E=210e9;poisson=0.2;h=0.001;F=10^6;
d_length=0.1;lengthx=2;lengthy=1;
Node=[];
for i=0:d_length:lengthx
 for j=0:d_length:lengthy
     node_num=node_num+1;
     Node=[Node;node_num,i,j];
 end
end
Num_eles=0;
element=[];
%该问题按直角三角形单元划分，单元坐标有两种形式，以下为左下角为直角的单元形式。
for i=d_length:d_length:lengthx
 for j=d_length:d_length:lengthy
 Num_eles=Num_eles+1;
 node_local=[i-d_length,j-d_length;i,j-d_length;i-d_length,j];
 node_list=[];
%已经有了单元的三个节点坐标，需要找到对应的节点编号
 for k=1:3
 row_x=find(abs(Node(:,2)-node_local(k,1))<1e-6);
 row_y=find(abs(Node(:,3)-node_local(k,2))<1e-6);
 num_node=intersect(row_x,row_y);
 %这样就找到了与 Node_info 中与节点横纵坐标均相等的节点编号了
 node_list=[node_list;num_node];
 %3x1 的矩阵，单元三个节点编号均找到
 end 
 element=[element;Num_eles,node_list'];
 end
end
nodes_num=size(Node,1);
eles_num=size(element,1);
K=zeros(2*nodes_num);
B=[];
for i=1:eles_num
 node_matrix=[element(i,2),Node(element(i,2),2),Node(element(i,2),3);
 element(i,3),Node(element(i,3),2),Node(element(i,3),3);
 element(i,4),Node(element(i,4),2),Node(element(i,4),3)];
 %单元节点矩阵，其第一列为节点编号，第二列为节点横坐标，第三列为节点纵坐标
 [Ke,matmtx,kinmtx2]=single_triangular(node_matrix,E,poisson,h); %求解单元刚度矩阵
 B=[B;kinmtx2]; 
 j=node_matrix(1,1);
 k=node_matrix(2,1);
 m=node_matrix(3,1);
 num=[2*j-1,2*j,2*k-1,2*k,2*m-1,2*m];
 for n1=1:6
 for n2=1:6
 K(num(n1),num(n2))=K(num(n1),num(n2))+Ke(n1,n2);  %刚度矩阵组装
 end
 end
end 
load_node_matrix=zeros(size(Node,1)*2,1);
row1_x=find(Node(:,2)==2);
row1_y=find(Node(:,3)==0);
num1=intersect(row1_x,row1_y);  %找到（2，0）点
row2_x=find(Node(:,2)==1);
row2_y=find(Node(:,3)==0);
num2=intersect(row2_x,row2_y);  %找到（1，0）点
row3_x=find(Node(:,2)==2);
row3_y=find(Node(:,3)==1);
num3=intersect(row3_x,row3_y);  %找到（2，1）点
row4_x=find(Node(:,2)==1);
row4_y=find(Node(:,3)==1);
num4=intersect(row4_x,row4_y);  %找到（1，1）点
load_node_matrix(2*num1)=-F/2;  %给四个节点赋值力
load_node_matrix(2*num3)=-F/2;
load_node_matrix(2*num2)=2*F;
load_node_matrix(2*num4)=-2*F;
num=find(Node(:,2)<1e-9);
KK=K;
ff=load_node_matrix;
for i=1:size(num,1)
 r=num(i);
 KK(2*r-1,:)=0;
 KK(:,2*r-1)=0;
 KK(2*r-1,2*r-1)=1;
 KK(2*r,:)=0;
 KK(:,2*r)=0;
 KK(2*r,2*r)=1;
 ff(2*r-1)=0;
 ff(2*r)=0;
end
d=pinv(KK)*ff;
d1=[];
d2=[];
for i=1:231
 d1=[d1;d(2*i-1)];
 d2=[d2;d(2*i)];
end
d11=d1;  %得到水平方向各节点位移
d22=d2;  %得到竖直方向各节点位移
X1=linspace(min(Node(:,2)),max(Node(:,2)),200);  %将横纵坐标均分成 200 份
Y1=linspace(min(Node(:,3)),max(Node(:,3)),200);
[X,Y,d1]=griddata(Node(:,2),Node(:,3),d1,X1',Y1,'natural');  %对二维或三维散点数据插值,插值方法为'natural'
[X,Y,d2]=griddata(Node(:,2),Node(:,3),d2,X1',Y1,'v4');       %插值方法为'v4'
subplot(4,1,1);
contourf(X,Y,d1);
colorbar  %分布柱图
title('水平方向位移分布等高图/m')
box on
subplot(4,1,2);
pcolor(X,Y,d1);
shading interp  %使色彩平滑过渡
colorbar  %分布柱图 
axis equal  %将横轴纵轴的定标系数设成相同值
title('水平方向位移分布云图/m')
box on;
x1_row=find(abs(Node(:,2)-2)<1e-5);  %找到右上角和右下角的节点标号
y1_row=find(abs(Node(:,3)-1)<1e-5);
y2_row=find(abs(Node(:,3))<1e-5);
node_num1=intersect(x1_row,y1_row);
node_num2=intersect(x1_row,y2_row);
x2_row=find(abs(Node(:,2)-0)<1e-5);  %找到左上角和左下角的节点标号
y3_row=find(abs(Node(:,3)-1)<1e-5);
y4_row=find(abs(Node(:,3))<1e-5);
node_num3=intersect(x2_row,y3_row);
node_num4=intersect(x2_row,y4_row);
text_x2y1=d11(node_num1);
text_x2y2=d11(node_num2);
text_x2y3=d22(node_num3);
text_x2y4=d22(node_num4);
text_1=['坐标','(',num2str(Node(node_num1,2)),',',num2str(Node(node_num1,3)),')',' 位移为',num2str(text_x2y1)];
text_2=['坐标','(',num2str(Node(node_num2,2)),',',num2str(Node(node_num2,3)),')',' 位移为',num2str(text_x2y2)];
text_3=['坐标','(',num2str(Node(node_num3,2)),',',num2str(Node(node_num3,3)),')',' 位移为',num2str(text_x2y3)];
text_4=['坐标','(',num2str(Node(node_num4,2)),',',num2str(Node(node_num4,3)),')',' 位移为',num2str(text_x2y4)];
text(Node(node_num1,2),Node(node_num1,3),text_1,'Color','magenta','FontSize',10);  %在云图中标注周边四个节点
text(Node(node_num2,2),Node(node_num2,3),text_2,'Color','magenta','FontSize',10);
text(Node(node_num3,2),Node(node_num3,3),text_3,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
text(Node(node_num4,2),Node(node_num4,3),text_4,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
subplot(4,1,3);
contourf(X,Y,d2);
colorbar  %分布柱图
title('竖直方向位移分布等高图/m')
box on;
subplot(4,1,4);
pcolor(X,Y,d2);
shading interp  %使色彩平滑过渡
colorbar  %分布柱图 
axis equal  %将横轴纵轴的定标系数设成相同值
title('竖直方向位移分布云图/m')
box on;
x1_row=find(abs(Node(:,2)-2)<1e-5);  %找到右上角和右下角的节点标号
y1_row=find(abs(Node(:,3)-1)<1e-5);
y2_row=find(abs(Node(:,3))<1e-5);
node_num1=intersect(x1_row,y1_row);
node_num2=intersect(x1_row,y2_row);
x2_row=find(abs(Node(:,2)-0)<1e-5);  %找到左上角和左下角的节点标号
y3_row=find(abs(Node(:,3)-1)<1e-5);
y4_row=find(abs(Node(:,3))<1e-5);
node_num3=intersect(x2_row,y3_row);
node_num4=intersect(x2_row,y4_row);
text_x2y1=d22(node_num1);
text_x2y2=d22(node_num2);
text_x2y3=d22(node_num3);
text_x2y4=d22(node_num4);
text_1=['坐标','(',num2str(Node(node_num1,2)),',',num2str(Node(node_num1,3)),')',' 位移为',num2str(text_x2y1)]; 
text_2=['坐标','(',num2str(Node(node_num2,2)),',',num2str(Node(node_num2,3)),')',' 位移为',num2str(text_x2y2)];
text_3=['坐标','(',num2str(Node(node_num3,2)),',',num2str(Node(node_num3,3)),')',' 位移为',num2str(text_x2y3)];
text_4=['坐标','(',num2str(Node(node_num4,2)),',',num2str(Node(node_num4,3)),')',' 位移为',num2str(text_x2y4)];
text(Node(node_num1,2),Node(node_num1,3),text_1,'Color','magenta','FontSize',10);  %在云图中标注周边四个节点
text(Node(node_num2,2),Node(node_num2,3),text_2,'Color','magenta','FontSize',10);
text(Node(node_num3,2),Node(node_num3,3),text_3,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
text(Node(node_num4,2),Node(node_num4,3),text_4,'Color','magenta','FontSize',10,'HorizontalAlignment','right');
sgtitle('实心结构位移分布图/m') ;