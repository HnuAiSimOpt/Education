%***********************************************************************
%                        一阶六面体单元主程序
%                         S230200242 刘齐耀 
%  demo3: 趴放，125N双集中力载荷, 网格模型信息: job_3.inp           
%***********************************************************************
function demo3()
%读取inp文件获得节点坐标信息Nodes及单元信息Elements
[Nodes, Elements] = Readmesh( 'job_3.inp' );
% 外力矩阵 Forces=[受力节点  受力方向(1,2,3分别代表x,y,z)  外力大小]  外力节点的编号在inp文件里面找
Forces=[1 2 125; 2 2 125]; 
%约束节点的编号在inp文件里面找
ConNumber=linspace(7201, 7272, 7272-7201+1);
%约束矩阵 Constraints=[强制位移节点  强制位移方向(1,2,3分别代表x,y,z)  强制位移大小]  
Constraints=zeros(size(ConNumber,2)*3,3);
for i=1:size(ConNumber,2)
Constraints(3*i-2:3*i,:)=[ConNumber(i) 1 0;ConNumber(i) 2 0;ConNumber(i) 3 0;];
end
E=2e5; %弹性模量
u=0.3;    %泊松比
%调用应变应力矩阵D
D=LinearIsotropicD(E,u);
U=StaticsSolver(E,u,Forces,Constraints,Nodes,Elements);
% 输出结果
OutputTXT = fopen('demo3_Results.txt','w'); %打开一个可写文件，用于写入计算结果
OutputResults(OutputTXT,Nodes,Elements,D,U)%调用输出结果文件
fclose(OutputTXT);
edit('demo3_Results.txt')
end