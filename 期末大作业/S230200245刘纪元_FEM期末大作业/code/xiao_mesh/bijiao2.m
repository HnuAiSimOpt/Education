clc;
clear
%导入有限元软件中的位移矩阵
uus = load('C:\Users\alein\Desktop\S200200239_FEM\code\xiao_mesh\txt\uu.txt');      
uu=[];
[m n]=size(uus);
for i=1:m
    uu(2*i-1, 1)=uus(i, 2);
    uu(2*i,1)=uus(i,3);
end
%导入matlab计算得到的位移矩阵
uu_fem= load('C:\Users\alein\Desktop\S200200239_FEM\code\xiao_mesh\txt\U_fem.txt');
%导入有限元软件中节点应力矩阵
s_node= load('C:\Users\alein\Desktop\S200200239_FEM\code\xiao_mesh\txt\s_node.txt');
%导入matlab计算得到的节点应力矩阵
stress_node= load('C:\Users\alein\Desktop\S200200239_FEM\code\xiao_mesh\txt\stress_node.txt');
%位移计算结果与软件分析结果的差值
difference_u=uu_fem-uu;
uu1=uu;
%x轴应力两种结果的差值
difference_sx=stress_node-s_node(:,2:4);
stress_x=stress_node(:,1);
s_x=s_node(:,2);
%合位移误差矩阵
error_u=[];
%x轴应力误差矩阵
error_s=[];
%软件分析结果与matlab计算结果的比值
e_u=[];
e_s=[];
error_u_abs_r=zeros(2*m,1);
for j=1:2*m
    error_u(j,1)=difference_u(j,1)/uu_fem(j,1);
    e_u(j,1)=uu1(j,1)/uu_fem(j,1);
    error_u_abs=abs(error_u);
   % if error_u_abs(j)==1
    %   error_u_abs_r(j)=0;
   % else
    %    error_u_ab_r(j)= error_u_abs(j);
    %end
end
for i=1:m
    error_s(i,1)=difference_sx(i,1)/stress_x(i,1);
    e_s(i,1)=s_x(i,1)/stress_x(i,1);
    error_s_abs=abs(error_s);
end
%计算结果与软件分析结果误差的均值
a=0;
for i=1:2*m
    a=a+error_u_abs(i);
end
%mean_u=mean(error_u_abs);
mean_u=a/(2*m);
mean_s=mean(error_s_abs);