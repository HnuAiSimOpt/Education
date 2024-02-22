
function PostContour(Nodes,Elements,U,Component_val)
NodeNum = size(Nodes,1) ;  %�ڵ����
ElementNum = size(Elements,1) ;    %��Ԫ����
ElementNodeNum=4;  %ÿ����Ԫ�ڵ���
% �����ʼ����X Y Z������꣺value���ֵ,��ÿ����Ԫ���սڵ�������λ�����ͼ
X = zeros(ElementNodeNum,ElementNum) ;%4*��Ԫ����������Ԫ���ڵ�ĺ�����
Y = zeros(ElementNodeNum,ElementNum) ;
Z = zeros(ElementNodeNum,ElementNum) ;
value = zeros(ElementNodeNum,ElementNum) ;
%�жϾ������ͣ�λ�ƣ�Ӧ����Ӧ�䣩
if size(Component_val,1)>1%λ��
    for i=1:ElementNum
        nd=Elements(i,:);
        value(:,i) = Component_val(nd) ;
    end
else %Ӧ��Ӧ��
    %�Ƚ���ĥƽ��avg ���ֵ��ֵ���ڵ㴦��ȡƽ����
    Difference=max(Component_val)-min(Component_val);%ȫ���ϵ����ֵ-��Сֵ
    AVG=0.75;  % Ĭ����ֵ75%
    for i=1:1:NodeNum %�����ڵ㣬����Ӧ��ĥƽ
        TElements=Elements';%ת��Elements
        itemp=(TElements==i);%�����߼��жϣ�itemp:Ԫ��Ϊ0 1�ľ���
        Cut=max(Component_val(1,itemp))-min(Component_val(1,itemp));%�ýڵ��Ӧ����Ӧ�䣩��
        if 0<Cut&&Cut<=AVG*Difference(1)%�ж��Ƿ�������ֵ����
            Component_val(1,itemp)=mean(Component_val(1,itemp));%����Ӧ����Ӧ�䣩ĥƽ
        end
    end
    value=reshape(Component_val,ElementNodeNum,ElementNum);%��Component��ֵ����value
end

% ���Ʊ��κ���ͼ
newNodes=Nodes';
newNodes=newNodes(:);
DeformationCoefficient=5.0e2; %���ηŴ�ϵ��
newNodes=newNodes+DeformationCoefficient*U;
newNodes=reshape(newNodes,[3,size(Nodes,1)]);
newNodes=newNodes';
% ���嵥Ԫÿ����Ľڵ�˳��(˳ʱ�����ʱ��)
if ElementNodeNum == 4  % C3D4��Ԫ
    fm = [1 2 4;2 3 4;3 1 4;1 3 2]; 
elseif ElementNodeNum==6 %C3D6��Ԫ
    fm=[1 2 3 1;4 5 6 4;1 2 5 4;1 3 6 4;2 3 6 5];
elseif ElementNodeNum == 8  % C3D8��Ԫ
    fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
elseif ElementNodeNum == 20  % C3D20��Ԫ
    fm = [1,9,2,10,3,11,4,12;5,13,6,14,7,15,8,16;1,9,2,18,6,13,5,17;
        2,10,3,19,7,14,6,18;3,11,4,20,8,15,7,19;1,17,5,16,8,20,4,12]; 
end
xyz = cell(1,ElementNum) ;
profile = xyz ;
%���ڵ�λ�ø������Ӧֵ
for e=1:ElementNum %ѭ����ȡÿ����Ԫ�ڽڵ������
    nd=Elements(e,:);
    X = newNodes(nd,1) ;
    Y = newNodes(nd,2) ;
    Z = newNodes(nd,3) ;
    xyz{e} = [X  Y Z] ;
    profile{e} = value(:,e);
end
% ����patch�������л�ͼ
figure
cellfun(@patch,repmat({'Vertices'},1,ElementNum),xyz,.......%�ö������Ƭ�ķ�ʽ��ʾ������
    repmat({'Faces'},1,ElementNum),repmat({fm},1,ElementNum),......
    repmat({'FaceVertexCdata'},1,ElementNum),profile,......
    repmat({'FaceColor'},1,ElementNum),repmat({'interp'},1,ElementNum));
view(3);
rotate3d on;
axis off; %����ʾ������
colormap(jet);
caxis([min(Component_val),max(Component_val)]);
t1=caxis;
t1=linspace(t1(1),t1(2),13);
colorbar('ytick',t1,'Location','westoutside');
axis equal;
end