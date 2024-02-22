
%%---------------------------comment-----------------------------------%%
% 姓名 游婉婷                                                           
% 学号 S230200182                                                      
% purpose：采用线性三角形解决中空矩形薄板受集中力作用的应力应变问题，        
% 比较不挖空与中空位移各节点的区别并绘制水平和竖直方向上的位移云图及应力云图，
% 从位移等高图可看出明显差异。后面的sigma_x、sigma_y云图画的是不挖孔结构
% 具体参数设置及条件见readme.docx
%%----------------------------end--------------------------------------%%

%% 程序开始
clear all
% 定义结构参数
lengthx=2;       % 薄板长度
lengthy=1;       % 薄板宽度
d_length=0.1;    % 网格单元长度
h=0.001;         % 薄板厚度
E=210e9;         % 弹性模量
poisson=0.2;     % 泊松比
F=10^6;          % 施加外力值
% 三角形网格划分网格
[element,node]=mesh_grid(lengthx,lengthy,d_length); 
% 求解单元刚度矩阵并组装成总刚度矩阵
[matmtx,K,B]=assemble_K(node,element,E,poisson,h);
% 根据载荷分布将力赋给节点
load_node_matrix=load_to_node(node,F);
% 施加边界条件，消除总刚度矩阵的奇异性（确保位移有解） 
[KK,ff]=feaplyc2(node,K,load_node_matrix); 
% 采用违逆求解各单元节点的位移
d=pinv(KK)*ff;  
% 求解单元应力
[sigma_x,sigma_y,tao_xy]=element_stress(B,element,matmtx,d);  
% 绘制位移分布云图
plot_d(d,node);  
% 绘制 sigma_x 应力云图
plot_sigma_x(sigma_x,node,element);  
% 绘制 sigma_y 应力云图
plot_sigma_y(sigma_y,node,element); 

%% 程序结束
disp('程序运行结束，感谢亲爱的师兄检查程序！');