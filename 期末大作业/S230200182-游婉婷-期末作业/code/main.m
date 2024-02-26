clear all;
%定义参数
E=1.28e11;
poisson=0.3;
F=7.5e3;
load_k=8.0e4;  
kFileStr='mesh.k';   % 读取k文件
overview=importdata(kFileStr,',');
element_solid=Solve_element(overview);   % 求解单元信息
node=Solve_node(overview);   % 求解节点坐标
[K]=Solve_K(poisson,element_solid,node,E);   % 求解刚度矩阵
uni_load=Hexahedral3D8Node_Load(node,load_k,F);   % 得到载荷矩阵
[KK,ff]=Hexahedral3D8Node_Feaplyc2(node,K,uni_load);   % 保证刚度矩阵非奇异
[LL, UU]=lu(KK);
utemp=LL\ff;
d=UU\utemp;   % 求解节点位移
[sx,sy,sz,sxy,syz,szx]=Solve_stress(element_solid,node,E,poisson,d);   % 求解应力
Plot_d(node,element_solid,d)   % 绘制横截面变形示意图
Createplt_d(node,element_solid,d)   % 生成位移plt文件
Createplt_Sx(node,element_solid,sx)   % 生成x方向应力plt文件
Createplt_Sy(node,element_solid,sy)   % 生成y方向应力plt文件
Createplt_Sz(node,element_solid,sz)   % 生成z方向应力plt文件
Createplt_Sxy(node,element_solid,sxy)   % 生成xy平面切应力plt文件
Createplt_Syz(node,element_solid,syz)   % 生成yz平面切应力plt文件
Createplt_Szx(node,element_solid,szx)   % 生成zx平面切应力plt文件