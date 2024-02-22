%%姓名：吴靓  学号：S230200201
%程序可根据输入的任意四边形四个节点按给定离散单元数离散成三角形单元并计算节点位移、单元应变
%测试数据
%输入厚度:1
%输入弹性模量:1
%输入泊松比:0
%输入节点坐标      PPT例子：[0 0;0 -2;4 -1;4 0]
%输入x方向离散单元数 :16 
%输入y方向离散单元数 :8
%输入节点荷载      例[0;0;0;0;0;-1;0;0]
%给定边界条件      同载荷输入 例   [0;0;1;1;1;1;0;0]
clc;
clear ;
%%输入基本参数 
h=input("请输入厚度:");               %设置厚度
E=input("请输入弹性模量:");           %设置弹性模量
u=input("请输入泊松比:");             %设置泊松比
xy=input("请依次输入节点坐标x、y:");   %输入节点坐标      PPT例子：[0 0;0 -2;4 -1;4 0]
lx=input("请输入x方向离散单元数:");    %输入x方向离散单元数               
ly=input("请输入y方向离散单元数:");    %输入y方向离散单元数 
element=2*lx*ly;                       %划分的三角形单元数量/个
nn=(1+lx)*(ly+1);                    %节点数（number of nodes,nn)

    x(1)=xy(1,1);y(1)=xy(1,2);
    x(2)=xy(2,1);y(2)=xy(2,2);
    x(3)=xy(3,1);y(3)=xy(3,2);
    x(4)=xy(4,1);y(4)=xy(4,2);
    
%将节点重新排序(最低端节点为1，逆时针编号）
for i=1:3
    if y(i)>y(i+1)
        tx=x(i); ty=y(i);
        x(i)=x(i+1); y(i)=y(i+1); 
        x(i+1)=tx; y(i+1)=ty;
    end
    if x(i)<x(i+1)
        tx=x(i); ty=y(i); 
        x(i)=x(i+1); y(i)=y(i+1); 
        x(i+1)=tx; y(i+1)=ty; 
    end
end
%平移图像，将节点1变为原点（0，0）
for i=2:4
    x(i)=x(i)-x(1);y(i)=y(i)-y(1);
end
x(1)=0;y(1)=0;
x(5)=x(1);y(5)=y(1);      %便于计算边长
%%确定各单元的节点
%计算边长
length=zeros(4,1);    %存放边长顺序按节点1-2、2-3、3-4、4-1
for i=1:4
    length(i)=sqrt((x(i)-x(i+1))^2+(y(i)-y(i+1))^2);
end

xl=[];    %存放左侧边界4的节点
sinl=(y(4)-y(1))/length(4);
cosl=(x(4)-x(1))/length(4);
xr=[];    %存放右侧边界2的节点
sinr=(y(3)-y(2))/length(2);
cosr=(x(3)-x(2))/length(2);
for i=1:ly+1
    xl=[xl;x(1)+(i-1)*cosl*length(4)/ly y(1)+(i-1)*sinl*length(4)/ly];
    xr=[xr;x(2)+(i-1)*cosr*length(2)/ly y(2)+(i-1)*sinr*length(2)/ly];
end

x0=[];        %存放离散网格后的三角形单元节点坐标 
for i=1:ly+1
    for j=1:lx+1
        lengthi=sqrt((xl(i,1)-xr(i,1))^2+(xl(i,2)-xr(i,2))^2);
        sini=(xr(i,2)-xl(i,2))/lengthi;
        cosi=(xr(i,1)-xl(i,1))/lengthi;
        x0=[x0;xl(i,1)+(j-1)*cosi*lengthi/lx xl(i,2)+(j-1)*sini*lengthi/lx];
    end
end

%网格节点可视化
plot(x0(:,1),x0(:,2),'b.','MarkerSize',18)    
hold on

en=[];   %依次输入组成各单元的节点代号 
for i=1:ly
    for j=1:lx
        en=[en; (lx+1)*(i-1)+j (lx+1)*i+j (lx+1)*(i-1)+j+1;];
        en=[en; (lx+1)*i+j (lx+1)*i+j+1 (lx+1)*(i-1)+j+1;];
    end
end

%网格可视化
for i=1:ly   %竖直线、横线
    for j=1:lx
    line([x0((i-1)*(lx+1)+j,1),x0((i-1)*(lx+1)+j+1,1)],[x0((i-1)*(lx+1)+j,2),x0((i-1)*(lx+1)+j+1,2)],'linewidth',1)
    hold on
    line([x0((i-1)*(lx+1)+j,1),x0((i-1)*(lx+1)+j+lx+1,1)],[x0((i-1)*(lx+1)+j,2),x0((i-1)*(lx+1)+j+lx+1,2)],'linewidth',1)
    end
end
for i=1:ly   %斜线
    for j=1:lx
        line([x0((i-1)*(lx+1)+j+1,1),x0((i-1)*(lx+1)+j+lx+1,1)],[x0((i-1)*(lx+1)+j+1,2),x0((i-1)*(lx+1)+j+lx+1,2)],'linewidth',1)
        hold on
    end
end
for i=1:lx   %上边界
    line([x0(i+(1+lx)*ly,1),x0(1+i+(1+lx)*ly,1)],[x0(i+(1+lx)*ly,2),x0(1+i+(1+lx)*ly,2)],'linewidth',1)
    hold on
end
for i=1:ly   %右边界
    line([x0(i*(lx+1),1),x0((i+1)*(lx+1),1)],[x0(i*(lx+1),2),x0((i+1)*(lx+1),2)],'linewidth',1)
    hold on
end

%%计算刚度矩阵
K=zeros(2*nn,2*nn);        %定义初始刚度矩阵
p=[];                      %存放三角形单元的节点坐标
 for m=1:element
	for i=1:3
	for j=1:2
		p(i,j)=x0(en(m,i),j);
    end
    % z存放节点代号
    z(i)=en(m,i);
    end
%计算该单元的刚度矩阵

    x(1)=p(1,1);y(1)=p(1,2);
    x(2)=p(2,1);y(2)=p(2,2);
    x(3)=p(3,1);y(3)=p(3,2);

ZM(m,:)=z;   %存放三角形单元的节点代号，用于单元应力计算
   
     b1=y(2)-y(3);
     b2=y(3)-y(1);
     b3=y(1)-y(2);
     
     c1=x(3)-x(2);
     c2=x(1)-x(3);
     c3=x(2)-x(1);

     %计算三角形单元的面积Ae
     Ae=0.5*((x(2)*y(3)-x(3)*y(2))+(y(2)-y(3))*x(1)+(x(3)-x(2))*y(1));
     B=zeros(3,6);    %定义B矩阵
     B(1,1)=b1;       %对矩阵赋值
     B(1,3)=b2;
     B(1,5)=b3;

     B(2,2)=c1;
     B(2,4)=c2;
     B(2,6)=c3;

     B(3,1)=c1;
     B(3,2)=b1;
     B(3,3)=c2;
     B(3,4)=b2;
     B(3,5)=c3;
     B(3,6)=b3;

     B=B/(2*Ae);
%      BM=zeros(3,12);
     BM(:,6*m-5:6*m)=B;    %存放三角形单元的B矩阵，用于单元应力计算
     %计算D矩阵
     D=zeros(3);
     D(1,1)=1;
     D(1,2)=u;
     D(2,1)=u;
     D(2,2)=1;
     D(3,3)=(1-u)/2;
     D=E*D/(1-u*u);
     TB=B';
     ke=h*Ae*TB*D*B;
     %将三角形单元刚度矩阵ke赋值到整体刚度矩阵K中
     for i=1:3
         for j=1:3
     K(2*z(i)-1,2*z(j)-1)=K(2*z(i)-1,2*z(j)-1)+ke(2*i-1,2*j-1);
     K(2*z(i)-1,2*z(j))=K(2*z(i)-1,2*z(j))+ke(2*i-1,2*j);
     K(2*z(i),2*z(j)-1)=K(2*z(i),2*z(j)-1)+ke(2*i,2*j-1);
     K(2*z(i),2*z(j))=K(2*z(i),2*z(j))+ke(2*i,2*j);
         end
     end
 end
%  K=K*16;
 disp("整体刚度矩阵：")
 disp(K);

 %%输入等效荷载  (按照编号后的数顺序给定x,y方向荷载）
 f=zeros(2*nn,1);          
 finput=input("请依次输入各节点x、y方向荷载:");    %输入节点荷载 例[0;0;0;0;0;-1;0;0]
 for i=1:4
     if finput(2*i-1)~=0
         if 2*i-1==11
     elseif 2*i-1==3
             f(2*(lx+1)-1)=finput(2*i-1);
     elseif 2*i-1==5
             f(2*(lx+1)*(ly+1)-1)=finput(2*i-1);
     elseif 2*i-1==7
             f(2*((lx+1)*ly+1)-1)=finput(2*i-1);
         end
     end
             
 if finput(2*i)~=0
             if 2*i==2
             f(2)=finput(2*i);
             elseif 2*i==4
             f(2*(lx+1))=finput(2*i);
             elseif 2*i==6
             f(2*(lx+1)*(ly+1))=finput(2*i);
             elseif 2*i==8
             f(2*((lx+1)*ly+1))=finput(2*i);
             end
 end
 end
  
 %%给定边界条件
 U=ones(2*nn,1);   % 同载荷输入 例   [0;0;1;1;1;1;0;0]
 Uinput=input("请输入边界条件(位移为0时输入0，不为0时输入1):");
      if Uinput(1)==0&&Uinput(3)==0
         for j=1:lx+1
             U(2*j-1)=0;
         end
     elseif Uinput(3)==0&&Uinput(5)==0
         for j=1:ly+1
             U(2*(lx+1)*j-1)=0;
         end
     elseif Uinput(5)==0&&Uinput(7)==0
         for j=1:lx+1
             U(2*(ly*(lx+1)+j)-1)=0;
         end
     elseif Uinput(7)==0&&Uinput(1)==0
         for j=1:ly+1
             U(2*((lx+1)*(j-1)+1)-1)=0;
         end
      end

     if Uinput(2)==0&&Uinput(4)==0
         for j=1:lx+1
             U(2*j)=0;
         end
     elseif Uinput(4)==0&&Uinput(6)==0
         for j=1:ly+1
             U(2*(lx+1)*j)=0;
         end
     elseif Uinput(6)==0&&Uinput(8)==0
         for j=1:lx+1
             U(2*(ly*(lx+1)+j))=0;
         end
     elseif Uinput(8)==0&&Uinput(2)==0
         for j=1:ly+1
             U(2*((lx+1)*(j-1)+1))=0;
         end
     end
%将固定的节点相关刚度矩阵值变为0
 for i=1:2*nn
	if U(i)==0
 for m=1:2*nn
	K(m,i)=0;
	K(i,m)=0;
 end
	f(i)=0;
    K(i,i)=1;
   end
 end

 %%计算节点位移
 de=K\f;
 disp('节点位移')
 disp(de)

 %%计算各单元应变应力
     %获取该单元节点位移
     for i=1:3
         for m=1:element             
         dei(m,(2*i-1))=de(2*ZM(m,i)-1);
         dei(m,(2*i))=de(2*ZM(m,i));
         end
     end   
    
 
 %计算三角形单元应力
 %三角形单元应力
 for m=1:element
     stress(m,:)=D*BM(1:3,6*m-5:6*m)*dei(m,:)';
     strain(m,:)=stress(m,:)/E;
 end
 disp("应力计算结果：")
 for m=1:element
     disp('单元对应输入节点代号 ')
     disp(ZM(m,:))
     disp('单元应力')
     disp(stress(m,:)')
 end