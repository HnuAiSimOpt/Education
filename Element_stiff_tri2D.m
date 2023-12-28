function Ke=Element_stiff_tri2D(D,coords,h) % ,bodyforce)
% D -  弹性矩阵,形心处
% coords - 3个结点的坐标
% Ke - 单元刚度矩阵
% Fe -  体积力引起的单元节点力向量
% h - 单元厚度
% bodyforce -体积力，形心处

area=0.5*det([ones(3,1),coords]); %三角形面积
index=[1,2,3,1,2];
b=coords(index(2:4),2)-coords(index(3:5),2);
c=-coords(index(2:4),1)+coords(index(3:5),1);
% b=zeros(3,1);c=zeros(3,1);
% for i=1:3
%     b(i)=coords(index(i+1),2)-coords(index(i+2),2);
%     c(i)=-coords(index(i+1),1)+coords(index(i+2),1);
% end

B=zeros(3,6);
B(1,1:2:6)=b';B(2,2:2:6)=c';
B(3,1:2:6)=c';B(3,2:2:6)=b';
B=B/2/area;

Ke=B'*D*B*h*area;
% Fe=zeros(6,1);
% Fe(1:2:6)=h*area*bodyforce(1)*[1;1;1]/3;
% Fe(2:2:6)=h*area*bodyforce(2)*[1;1;1]/3;



