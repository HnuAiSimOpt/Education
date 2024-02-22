%%%输出结果程序(输出云图)
%  OutputTXT输出信息保存位置
%  Nodes节点坐标信息 
%  Elements单元信息
%  D计算各向同性线弹性材料应力-应变矩阵
%  U位移矩阵
function OutputResults(Nodes,Elements,D,U)
NodeCount = size(Nodes,1) ;  % 节点个数
ElementCount= size(Elements,1);%单元个数
ElementNodeCount=8;% 每个单元节点数
Dof=3; 
%计算高斯点、节点应力应变值
[NodeStrain,NodeStress,~,~]=CalculateStrainAndStress(U,D,Nodes,Elements);
%求得MISES应力矩阵
MISES=zeros(1,ElementCount*ElementNodeCount);
for I=1:ElementCount*ElementNodeCount
    MISES(I)=sqrt(0.5)*sqrt((NodeStress(1,I)-NodeStress(2,I))^2+(NodeStress(1,I)-NodeStress(3,I))^2+....
        (NodeStress(2,I)-NodeStress(3,I))^2+6*(NodeStress(4,I)^2+NodeStress(5,I)^2+NodeStress(6,I)^2));
end
%求得Umag位移矩阵
Umag=zeros(NodeCount,1);
for i=1:NodeCount
    Umag(i)=sqrt(U(3*i-2)^2+U(3*i-1)^2+U(3*i)^2);
end
%绘制变形前网格
  for i=1:1:size(Elements,1)
    points=Nodes(Elements(i,:),:);
    mesh=1:1:8;%网格信息
    %六面体单元节点坐标
    vertices_matrix = [points(mesh(1,:),1),points(mesh(1,:),2),points(mesh(1,:),3)];
    %六面体单元节点顺序
    faces_matrix= [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];%给出每个面的节点序号
    patch('vertices', vertices_matrix,'faces',faces_matrix,'facecolor','g');
    view(3);hold on%绘图
  end
axis equal
alpha(1);
 %输出位移、应力、应变云图
for i=1:Dof %输出U1-U3
    PlotContour(Nodes,Elements,U,U(i:3:size(U,1)))
    title(['U',num2str(i)]);%绘制云图标题
end
for i=1:size(NodeStress,1)
    %输出主应力、应变云图 S11-S33，E11-E33
    if i<4
        PlotContour(Nodes,Elements,U,NodeStrain(i,:))
        title(['线应变E',num2str(i),num2str(i)])
        PlotContour(Nodes,Elements,U,NodeStress(i,:))
        title(['正应力S',num2str(i),num2str(i)])
    %输出S12 S23 E12 E23云图
    elseif    i<6
        PlotContour(Nodes,Elements,U,NodeStrain(i,:))
        title(['角应变E',num2str(i-3),num2str(i-2)])
        PlotContour(Nodes,Elements,U,NodeStress(i,:))
        title(['切应力S',num2str(i-3),num2str(i-2)])
    %输出S13 E13云图
    else
        PlotContour(Nodes,Elements,U,NodeStrain(i,:))
        title(['角应变E',num2str(i-5),num2str(i-3)])
        PlotContour(Nodes,Elements,U,NodeStress(i,:))
        title(['切应力S',num2str(i-5),num2str(i-3)])
    end
end   
%绘制合位移云图
PlotContour(Nodes,Elements,U,Umag)
title('合位移Umag')
%绘制合应力云图
PlotContour(Nodes,Elements,U,MISES)
title('合应力MISES')
end