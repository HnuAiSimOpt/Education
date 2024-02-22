
function PostContour(Nodes,Elements,U,Component_val)
NodeNum = size(Nodes,1) ;  %节点个数
ElementNum = size(Elements,1) ;    %单元个数
ElementNodeNum=4;  %每个单元节点数
% 矩阵初始化，X Y Z点的坐标：value点的值,对每个单元按照节点序号依次绘制云图
X = zeros(ElementNodeNum,ElementNum) ;%4*单元个数；各单元各节点的横坐标
Y = zeros(ElementNodeNum,ElementNum) ;
Z = zeros(ElementNodeNum,ElementNum) ;
value = zeros(ElementNodeNum,ElementNum) ;
%判断矩阵类型（位移，应力，应变）
if size(Component_val,1)>1%位移
    for i=1:ElementNum
        nd=Elements(i,:);
        value(:,i) = Component_val(nd) ;
    end
else %应力应变
    %先进行磨平（avg 积分点插值到节点处再取平均）
    Difference=max(Component_val)-min(Component_val);%全域上的最大值-最小值
    AVG=0.75;  % 默认阈值75%
    for i=1:1:NodeNum %遍历节点，进行应力磨平
        TElements=Elements';%转置Elements
        itemp=(TElements==i);%进行逻辑判断，itemp:元素为0 1的矩阵
        Cut=max(Component_val(1,itemp))-min(Component_val(1,itemp));%该节点的应力（应变）差
        if 0<Cut&&Cut<=AVG*Difference(1)%判断是否满足阈值条件
            Component_val(1,itemp)=mean(Component_val(1,itemp));%进行应力（应变）磨平
        end
    end
    value=reshape(Component_val,ElementNodeNum,ElementNum);%将Component的值赋给value
end

% 绘制变形后云图
newNodes=Nodes';
newNodes=newNodes(:);
DeformationCoefficient=5.0e2; %变形放大系数
newNodes=newNodes+DeformationCoefficient*U;
newNodes=reshape(newNodes,[3,size(Nodes,1)]);
newNodes=newNodes';
% 定义单元每个面的节点顺序(顺时针或逆时针)
if ElementNodeNum == 4  % C3D4单元
    fm = [1 2 4;2 3 4;3 1 4;1 3 2]; 
elseif ElementNodeNum==6 %C3D6单元
    fm=[1 2 3 1;4 5 6 4;1 2 5 4;1 3 6 4;2 3 6 5];
elseif ElementNodeNum == 8  % C3D8单元
    fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
elseif ElementNodeNum == 20  % C3D20单元
    fm = [1,9,2,10,3,11,4,12;5,13,6,14,7,15,8,16;1,9,2,18,6,13,5,17;
        2,10,3,19,7,14,6,18;3,11,4,20,8,15,7,19;1,17,5,16,8,20,4,12]; 
end
xyz = cell(1,ElementNum) ;
profile = xyz ;
%将节点位置赋予其对应值
for e=1:ElementNum %循环获取每个单元内节点的坐标
    nd=Elements(e,:);
    X = newNodes(nd,1) ;
    Y = newNodes(nd,2) ;
    Z = newNodes(nd,3) ;
    xyz{e} = [X  Y Z] ;
    profile{e} = value(:,e);
end
% 采用patch函数进行绘图
figure
cellfun(@patch,repmat({'Vertices'},1,ElementNum),xyz,.......%用多边形面片的方式显示各个面
    repmat({'Faces'},1,ElementNum),repmat({fm},1,ElementNum),......
    repmat({'FaceVertexCdata'},1,ElementNum),profile,......
    repmat({'FaceColor'},1,ElementNum),repmat({'interp'},1,ElementNum));
view(3);
rotate3d on;
axis off; %不显示坐标轴
colormap(jet);
caxis([min(Component_val),max(Component_val)]);
t1=caxis;
t1=linspace(t1(1),t1(2),13);
colorbar('ytick',t1,'Location','westoutside');
axis equal;
end