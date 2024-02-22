clear all;
first_time=cputime;
format short e % 设定输出类型
fprintf=fopen('input.txt ','rt'); % 打开输入数据文件，读入参数数据
nelement=fscanf(fprintf,'%d',1);% 单元个数
npiont=fscanf(fprintf,'%d',1);% 结点个数
nbccondit=fscanf(fprintf,'%d',1) % 受约束边界点数
nforce=fscanf(fprintf,'%d',1);% 结点荷载个数
young=fscanf(fprintf,'%e',1);% 弹性模量
poission=fscanf(fprintf,'%f',1);% 泊松比
thickness=fscanf(fprintf,'%f',1);% 厚度
nodes=fscanf(fprintf,'%d',[3,nelement])';% 单元定义数组（单元结点号）
ncoordinates=fscanf(fprintf,'%f',[2,npiont])'; % 结点坐标数组
force=fscanf(fprintf,'%f',[3,nforce])'; % 结点力数组（受力结点编号 , x 方向 ,y 方向）
fc=fopen('constraint.txt ','rt');
constraint=fscanf(fc,'%d',[3,nbccondit])'; % 约束信息（约束点， x 约束， y 约束）%有约束为 1，无约束为 0
kk=zeros(2*npiont,2*npiont); % 生成特定大小总体刚度矩阵并置 0
for i=1:nelement
    D= [1 poission 0; 
 poission 1 0; 
 0 0 (1-poission)/2]*young/(1-poission^2) %生成弹性矩阵 D
A=det([1 ncoordinates(nodes(i,1),1) ncoordinates(nodes(i,1),2); 
 1 ncoordinates(nodes(i,2),1) ncoordinates(nodes(i,2),2); 
 1 ncoordinates(nodes(i,3),1) ncoordinates(nodes(i,3),2)])/2 %计算当前单元的面积A
for j=0:2 
b(j+1)=ncoordinates(nodes(i,(rem((j+1),3))+1) ,2)-ncoordinates(nodes(i,(rem((j+2),3))+1),2); 
c(j+1)=-ncoordinates(nodes(i,(rem((j+1),3))+1),1)+ncoordinates(nodes(i,(rem((j+2),3))+1),1);
end 
 B=[b(1) 0 b(2) 0 b(3) 0;
 0 c(1) 0 c(2) 0 c(3);
 c(1) b(1) c(2) b(2) c(3) b(3)]/(2*A); %生成应变矩阵 B
B1( :,:,i)=B;
S=D*B;%求应力矩阵S
nk=B'*S*thickness*A; % 求解单元刚度矩阵
 a=nodes(i,:); % 临时向量,用来记录当前单元的节点编号
 for j=1:3 
 for k=1:3 
 kk((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2)=kk((a(j)*2-1):a(j)*2,(a(k)*2-1):a(k)*2)+nk(j*2-1:j*2,k*2-1:k*2); 
 % 根据节点编号对应关系将单元刚度分块叠加到总刚度矩阵中
 end 
 end 
 end
%***************************************************************************** 
%将约束信息加入总体刚度矩阵（对角元素改一法）
 for i=1:nbccondit
 if constraint(i,2)==1 
 kk(:,(constraint(i,1)*2-1))=0; % 一列为零
 kk((constraint(i,1)*2-1),:)=0; % 一行为零
 kk((constraint(i,1)*2-1),(constraint(i,1)*2-1))=1; % 对角元素为 1 
 end
if constraint(i,3)==1 
 kk( :,constraint(i,1)*2)=0; % 一列为零
 kk(constraint(i,1)*2,:)=0; % 一行为零
 kk(constraint(i,1)*2 ,constraint(i,1)*2)=1; % 对角元素为 1 
 end 
 end 
%*****************************************************************************
%生成荷载向量
 loadvector(1:2*npiont)=0; % 总体荷载向量置零
 for i=1:nforce 
 loadvector((force(i,1)*2-1):force(i,1)*2)=force(i,2:3); 
 end
 %***************************************************************************** 
%求解内力
 displancement=kk\loadvector' % 计算节点位移向量
 edisplancement(1:6)=0; % 当前单元节点位移向量
 for i=1:nelement
 for j=1:3 
 edisplancement(j*2-1:j*2)=displancement(nodes(i,j)*2-1:nodes(i,j)*2); 
 % 取出当前单元的节点位移向量
 end 
 i ;
 stress=D*B1(:, :, i)*edisplancement'; % 求内力
 stress1(i,:)=[stress'];
 stress_x(i)=stress1(i,1);
 stress_y(i)=stress1(i,2);
 stress_xy(i)=stress1(i,3);
 sigma1(i)=0.5*(stress_x(i)+stress_y(i))+sqrt((stress_xy(i))^2+(0.5*(stress_x(i)-stress_y(i)))^2);
 sigma2(i)=0.5*(stress_x(i)+stress_y(i))-sqrt((stress_xy(i))^2+(0.5*(stress_x(i)-stress_y(i)))^2);
 stress_vonmises(i)=sqrt(0.5*((sigma1(i)-sigma2(i))^2+(sigma1(i)-0)^2+(0-sigma2(i))^2));
 dlmwrite('stress_vonmises.txt',stress_vonmises);
 dlmwrite('stress.txt',stress1);
 dlmwrite('displancement.txt',displancement);
 end
set(0,'defaultfigurecolor','w')
%画vonmiss应力云图
s1=max(stress_vonmises);%求出最大应力
s2=min(stress_vonmises);%求出最小应力
a=(s1-s2)/9;%将应力分成9份
stress_range=zeros(1,10);
stress_range(1)=s1;%stress_range(1)为最大应力
for i=2:10
    stress_range(i)=stress_range(i-1)-a;
end
range=stress_range;
figure(1);
color=jet(9);%将彩虹色分成9份
for i=1:size(nodes)
    ElementCoodinate=[ncoordinates(nodes(i,1),:)
                      ncoordinates(nodes(i,2),:)
                      ncoordinates(nodes(i,3),:)];
    x=ElementCoodinate(:,1);
    y=ElementCoodinate(:,2);
    s=stress_vonmises(i);
    %将应力大小与颜色对应
    if (range(1)>=s)&&(s>range(2))
        ColorSpec=color(9,:);
    elseif (range(2)>=s)&&(s>range(3))
        ColorSpec=color(8,:);
    elseif (range(3)>=s)&&(s>range(4))
        ColorSpec=color(7,:);
    elseif (range(4)>=s)&&(s>range(5))
        ColorSpec=color(6,:);
    elseif (range(5)>=s)&&(s>range(6))
        ColorSpec=color(5,:);
    elseif (range(6)>=s)&&(s>range(7))
        ColorSpec=color(4,:);
    elseif (range(7)>=s)&&(s>range(8))
        ColorSpec=color(3,:);
    elseif (range(8)>=s)&&(s>range(9))
        ColorSpec=color(2,:);
    else 
        ColorSpec=color(1,:);
    end
    fill(x,y,ColorSpec);%画出单元应力对应的颜色
    hold on          
end 
% 画应力云图标签
range=sort(range);
colormap(color);
c=colorbar;
c.TickLabels=(range);
c.Ticks=[0,1/9,2/9,3/9,4/9,5/9,6/9,7/9,8/9,9/9];
axis equal;
title('vonmiss应力云图');

%画sigmay应力云图
s1=max(stress_y);%求出最大应力
s2=min(stress_y);%求出最小应力
a=(s1-s2)/9;%将应力分成9份
stress_range=zeros(1,10);
stress_range(1)=s1;%stress_range(1)为最大应力
for i=2:10
    stress_range(i)=stress_range(i-1)-a;
end
range=stress_range;
figure(2);
color=jet(9);%将彩虹色分成9份
for i=1:size(nodes)
    ElementCoodinate=[ncoordinates(nodes(i,1),:)
                      ncoordinates(nodes(i,2),:)
                      ncoordinates(nodes(i,3),:)];
    x=ElementCoodinate(:,1);
    y=ElementCoodinate(:,2);
    s=stress_y(i);
    %将应力大小与颜色对应
    if (range(1)>=s)&&(s>range(2))
        ColorSpec=color(9,:);
    elseif (range(2)>=s)&&(s>range(3))
        ColorSpec=color(8,:);
    elseif (range(3)>=s)&&(s>range(4))
        ColorSpec=color(7,:);
    elseif (range(4)>=s)&&(s>range(5))
        ColorSpec=color(6,:);
    elseif (range(5)>=s)&&(s>range(6))
        ColorSpec=color(5,:);
    elseif (range(6)>=s)&&(s>range(7))
        ColorSpec=color(4,:);
    elseif (range(7)>=s)&&(s>range(8))
        ColorSpec=color(3,:);
    elseif (range(8)>=s)&&(s>range(9))
        ColorSpec=color(2,:);
    else 
        ColorSpec=color(1,:);
    end
    fill(x,y,ColorSpec);%画出单元应力对应的颜色
    hold on          
end 
% 画应力云图标签
range=sort(range);
colormap(color);
c=colorbar;
c.TickLabels=(range);
c.Ticks=[0,1/9,2/9,3/9,4/9,5/9,6/9,7/9,8/9,9/9];
axis equal;
title('sigma_y应力云图');

%画sigmax应力云图
s1=max(stress_x);%求出最大应力
s2=min(stress_x);%求出最小应力
a=(s1-s2)/9;%将应力分成9份
stress_range=zeros(1,10);
stress_range(1)=s1;%stress_range(1)为最大应力
for i=2:10
    stress_range(i)=stress_range(i-1)-a;
end
range=stress_range;
figure(3);
color=jet(9);%将彩虹色分成9份
for i=1:size(nodes)
    ElementCoodinate=[ncoordinates(nodes(i,1),:)
                      ncoordinates(nodes(i,2),:)
                      ncoordinates(nodes(i,3),:)];
    x=ElementCoodinate(:,1);
    y=ElementCoodinate(:,2);
    s=stress_x(i);
    %将应力大小与颜色对应
    if (range(1)>=s)&&(s>range(2))
        ColorSpec=color(9,:);
    elseif (range(2)>=s)&&(s>range(3))
        ColorSpec=color(8,:);
    elseif (range(3)>=s)&&(s>range(4))
        ColorSpec=color(7,:);
    elseif (range(4)>=s)&&(s>range(5))
        ColorSpec=color(6,:);
    elseif (range(5)>=s)&&(s>range(6))
        ColorSpec=color(5,:);
    elseif (range(6)>=s)&&(s>range(7))
        ColorSpec=color(4,:);
    elseif (range(7)>=s)&&(s>range(8))
        ColorSpec=color(3,:);
    elseif (range(8)>=s)&&(s>range(9))
        ColorSpec=color(2,:);
    else 
        ColorSpec=color(1,:);
    end
    fill(x,y,ColorSpec);%画出单元应力对应的颜色
    hold on          
end 
% 画应力云图标签
range=sort(range);
colormap(color);
c=colorbar;
c.TickLabels=(range);
c.Ticks=[0,1/9,2/9,3/9,4/9,5/9,6/9,7/9,8/9,9/9];
axis equal;
title('sigma_x应力云图');

%画sigmaxy应力云图
s1=max(stress_xy);%求出最大应力
s2=min(stress_xy);%求出最小应力
a=(s1-s2)/9;%将应力分成9份
stress_range=zeros(1,10);
stress_range(1)=s1;%stress_range(1)为最大应力
for i=2:10
    stress_range(i)=stress_range(i-1)-a;
end
range=stress_range;
figure(4);
color=jet(9);%将彩虹色分成9份
for i=1:size(nodes)
    ElementCoodinate=[ncoordinates(nodes(i,1),:)
                      ncoordinates(nodes(i,2),:)
                      ncoordinates(nodes(i,3),:)];
    x=ElementCoodinate(:,1);
    y=ElementCoodinate(:,2);
    s=stress_xy(i);
    %将应力大小与颜色对应
    if (range(1)>=s)&&(s>range(2))
        ColorSpec=color(9,:);
    elseif (range(2)>=s)&&(s>range(3))
        ColorSpec=color(8,:);
    elseif (range(3)>=s)&&(s>range(4))
        ColorSpec=color(7,:);
    elseif (range(4)>=s)&&(s>range(5))
        ColorSpec=color(6,:);
    elseif (range(5)>=s)&&(s>range(6))
        ColorSpec=color(5,:);
    elseif (range(6)>=s)&&(s>range(7))
        ColorSpec=color(4,:);
    elseif (range(7)>=s)&&(s>range(8))
        ColorSpec=color(3,:);
    elseif (range(8)>=s)&&(s>range(9))
        ColorSpec=color(2,:);
    else 
        ColorSpec=color(1,:);
    end
    fill(x,y,ColorSpec);%画出单元应力对应的颜色
    hold on          
end 
% 画应力云图标签
range=sort(range);
colormap(color);
c=colorbar;
c.TickLabels=(range);
c.Ticks=[0,1/9,2/9,3/9,4/9,5/9,6/9,7/9,8/9,9/9];
axis equal;
title('sigmaxy应力云图');

fclose(fc); % 关闭数据文件
fclose(fprintf); % 关闭数据文件


