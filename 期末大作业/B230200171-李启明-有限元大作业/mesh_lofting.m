function [Node,Element] = mesh_lofting(p0,temp_X,temp_Y,z)
%拉伸成三维模型
%p0：几何模型原点
%x：x方向生成的节点坐标
%y：y方向生成的节点坐标
%z：z方向生成的节点坐标

nx = size(temp_X,2)-1;%x方向划分网格数
ny =size(temp_Y,2)-1;%y方向划分网格数
nz = length(z)-1;%z方向划分网格数

temp_X = temp_X+p0(1);
temp_Y = temp_Y+p0(2);
z = z+p0(3);

%第一层[x-y]平面
zi = z(1);
[X,Y] = meshgrid(temp_X(1,:),temp_Y(1,:));
X = X';
Y = Y';
Zi = repmat(zi,[length(X(:)),1]);
Node = [X(:),Y(:),Zi];

num = size(Node,1);%每层的节点数
NumE = nx*ny; %每层的单元数
%总的单元数
NumElement = nx*ny*nz;
Element = zeros(NumElement,8);%单元连接节点[1 2 3 4 5 6 7 8]

count = 1;
for k = [2:length(z)]
    [X,Y] = meshgrid(temp_X(k,:),temp_Y(k,:));
    X = X';
    Y = Y';
    Zi = repmat(z(k),[length(X(:)),1]);%第k层的z坐标
    len = size(Node,1);
    %第k-1层的节点序列
    series1 = [len-num+1:len];   
    %每层循环找出单元连接节点
    for n = [1:NumE]
        %每一层中的相对单元号
        count_p = count - (k-2)*NumE;
         
        %第k个单元的第1,2,3,4,5,6,7,8个节点
         if mod(count_p,nx)==0
            num1 = series1((floor(count_p/nx)-1)*(nx+1)+nx);
        else
            num1 =series1(floor(count_p/nx)*(nx+1)+mod(count_p,nx));
        end
        num2 = num1+1;
        num3 = num2+nx+1;
        num4 = num3-1;
        %第k层
        num5 =num1+num;
        num6 = num5+1;
        num7 = num6+nx+1;
        num8 = num7-1;
        
        Element(count,1:8) = [num1,num2,num3,num4,num5,num6,num7,num8];
      
        count = count+1;
    end
    
    Node = [Node;[X(:),Y(:),Zi]];
end

for ie = [1:NumElement]
      num1 = Element(ie,1);
      num2 = Element(ie,2);
      num4 = Element(ie,4);
      num5 = Element(ie,5);
end

end














