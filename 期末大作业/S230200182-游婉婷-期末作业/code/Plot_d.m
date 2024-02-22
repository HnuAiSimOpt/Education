function []=Plot_d(node,element_solid,d)
eles_num=size(element_solid,1);
nodes_num=size(node,1);
A=find(abs(node(:,4)-60)<1e-6);   % 定义横截面点
b=0;
for i=1:eles_num
    if (ismember(element_solid(i,2),A)==1||ismember(element_solid(i,3),A)==1||ismember(element_solid(i,4),A)==1||ismember(element_solid(i,5),A)==1)
        b=b+1;
        I(b)=i;    % 找到横截面点对应的单元编号
    end
end
X=d(1:3:8928);Y=d(2:3:8928);    % 得到横截面上点的x,y方向位移
for i=1:nodes_num
    X(i)=node(i,2)+X(i)*8e5;   % 由于变形太小，将位移量乘以放大系数达到同一量级
end
for i=1:nodes_num
    Y(i)=node(i,3)+Y(i)*8e5;
end
figure('name','Fig','Position',[400,200,600,350])
for i=1:124
    patch(node(element_solid(I(i),2:5),2),node(element_solid(I(i),2:5),3),node(element_solid(I(i),2:5),4),'w','FaceColor','none','LineStyle','--','EdgeColor','b')   % 未变形网格用蓝色虚线表示
    hold on
    patch(X(element_solid(I(i),2:5)),Y(element_solid(I(i),2:5)),node(element_solid(I(i),2:5),4),'w','FaceColor','none','EdgeColor','r')   % 变形后网格用红色实线表示
end
legend('变形前', '变形后')