<<<<<<< HEAD

%本程序采用的是CPS4单元即四节点四边形单元，计算格式采用的是显式格式，计算方法采用全积分，
%即在每个时间步都利用控制方程计算出所有节点的加速度，进而得到每个时间步节点的速度、位移和应力。
clc;
%前处理程序，直接建立模型，得到节点信息，单元信息等
lengthx=4;%模型大小
lengthy=2;
emodule=5e+6;%弹性模量
possion=0.25;%泊松比
density=0.0005;%密度
c=0.0;%阻尼
numelex=50;%x、y方向单元数
numeley=25;
numele=numelex*numeley;
numnode=(numelex+1)*(numeley+1);
nglx=2;%x、y方向高斯积分点数
ngly=2;
nnel=4;
coornode=[zeros(2,numnode)];%节点坐标
for i=1:numelex+1
    for j=1:numeley+1
        coornode(1,(j-1)*(numelex+1)+i)=(i-1)*lengthx/numelex;
        coornode(2,(j-1)*(numelex+1)+i)=(j-1)*lengthy/numeley;
    end
end

%单元对应节点坐标编号
nodes=[zeros(numele,4)];
mat1=1:numnode;
mat2=reshape(mat1,numelex+1,numeley+1);
for i=1:numele
    h=ceil(i/numelex);
    l=rem(i,numelex);
    if(l==0)
        l=numelex;
    end
    mat3=mat2(l:l+1,h:h+1);
    mat4=mat3(:);
    nodes(i,1:4)=[mat4([4,3,1,2])'];
end

%初始化节点内外力，施加力边界条件
exforce=[zeros(2,numnode)];
 %在最右端所有节点施加y方向大小为-1的集中力 
 for i=1:numeley+1
     exforce(2,(numelex+1)*i)=-1;
 end
inforce=[zeros(2,numnode)];
stress=[zeros(3,numnode)];

%计算每个节点质量，利用的一致质量矩阵原理，即每个单元的质量均分给它的节点
mass1=[zeros(1,numnode)];
for e=1:numele
    x1=coornode(1,nodes(e,1));
    x2=coornode(1,nodes(e,2));
    x3=coornode(1,nodes(e,3));
    x4=coornode(1,nodes(e,4));
    y1=coornode(2,nodes(e,1));
    y2=coornode(2,nodes(e,2));
    y3=coornode(2,nodes(e,3));
    y4=coornode(2,nodes(e,4));
    A=abs(0.5*(x2*y3+x3*y1+x1*y2-x3*y2-x1*y3-x2*y1));%单元面积
    A=A+abs(0.5*(x3*y4+x4*y2+x2*y3-x4*y3-x2*y4-x3*y2));
    m=density*A/4;
    mass1(1,nodes(e,1))=mass1(1,nodes(e,1))+m;
    mass1(1,nodes(e,2))=mass1(1,nodes(e,2))+m;
    mass1(1,nodes(e,3))=mass1(1,nodes(e,3))+m;
    mass1(1,nodes(e,4))=mass1(1,nodes(e,4))+m;
end
mass=[mass1;mass1];%节点质量矩阵

%初始化节点的位移、速度、加速度。
disp=[zeros(2,numnode)];
vel=[zeros(2,numnode)];
a=[zeros(2,numnode)];

%设置时间步
numstep=1000;
time=[zeros(1,numstep)];
D=(emodule/(1-possion^2))*[1,possion,0;possion,1,0;0,0,(1-possion)/2];%材料矩阵（平面应力类型）
[point2,weight2]=glqd2(nglx,ngly);%得到高斯积分点位置矩阵和权系数矩阵
stab=zeros(4,4);
nnd=zeros(1,numnode);%每个节点被单元共享次数
  for i=1:numele
      for j=1:4
          nnd(nodes(i,j))=nnd(nodes(i,j))+1;
      end
  end
  numnd=[nnd;nnd;nnd];
  nxdisp=zeros(3,numstep);
  nydisp=zeros(3,numstep);
 
  %开始时间步循环
for step=1:numstep
    inforce=[zeros(2,numnode)];
    stress=[zeros(3,numnode)];
    mindelt=0.0001;%初始设置单个时间步持续时间
    for e=1:numele %每个时间步单元循环
        glstress=[zeros(4,3)];
        X=zeros(1,4);
        Y=zeros(1,4);
        for i=1:4
         X(i)=coornode(1,nodes(e,i))+disp(1,nodes(e,i));
         Y(i)=coornode(2,nodes(e,i))+disp(2,nodes(e,i));
        end
         [a1,b1,c1,d1,angleA,angleB,angleC,angleD]=angle(X,Y);%得到每个单元内角
         %判断单元是否失效
       if angleA >=pi || angleB >=pi ||angleC >=pi || angleD >=pi 
          temp='quadrilateral element is concave - FEM can not solve';
          break;
       end
       
       %根据每个单元的最小边长计算最小时间步持续时间
       delt=min([a1,b1,c1,d1]);
       delt=0.7*delt/sqrt(emodule/density);
        if (delt<mindelt)
            mindelt=delt;
        end
        
        pk=0;
        for  inty=1:ngly %积分点循环
                y=point2(inty,2);      
                wty=weight2(inty,2);     
            for intx=1:nglx
             x=point2(intx,1);         
            wtx=weight2(intx,1);  
                [shape,dhdr,dhds]=glsb(x,y);%得到自然坐标下每个积分点的形函数及其偏导数
                jacob2=jacob(nnel,dhdr,dhds,X,Y); %得到高斯点的雅各比矩阵
                detjacob=det(jacob2);                        
                invjacob=inv(jacob2);
                [dhdx,dhdy]=glb(nnel,dhdr,dhds,invjacob);%得到物理坐标下节点形函数及其偏导数
                glB=feglB(nnel,dhdx,dhdy);%组装物理坐标下节点的B矩阵（应变矩阵）
                q=[zeros(8,1)];
                q(1,1)=disp(1,nodes(e,1));
                q(3,1)=disp(1,nodes(e,2));
                q(5,1)=disp(1,nodes(e,3));
                q(7,1)=disp(1,nodes(e,4));
                q(2,1)=disp(2,nodes(e,1));
                q(4,1)=disp(2,nodes(e,2));
                q(6,1)=disp(2,nodes(e,3));
                q(8,1)=disp(2,nodes(e,4));
                glnstrain=glB*q;%得到自然坐标的节点应变
                glnstress=D*glnstrain;%得到自然坐标的节点应力
                 pk=pk+1;
                 stab(pk,:)=shape;%组装应力转换矩阵
                glstress(pk,:)=glnstress';
                ninforce=detjacob*glB'*glnstress*wtx*wty;%得到物理坐标下节点的内力
                inforce(1,nodes(e,1))=inforce(1,nodes(e,1))+ninforce(1,1);
                inforce(2,nodes(e,1))=inforce(2,nodes(e,1))+ninforce(2,1);
                inforce(1,nodes(e,2))=inforce(1,nodes(e,2))+ninforce(3,1);
                inforce(2,nodes(e,2))=inforce(2,nodes(e,2))+ninforce(4,1);
                inforce(1,nodes(e,3))=inforce(1,nodes(e,3))+ninforce(5,1);
                inforce(2,nodes(e,3))=inforce(2,nodes(e,3))+ninforce(6,1);
                inforce(1,nodes(e,4))=inforce(1,nodes(e,4))+ninforce(7,1);
                inforce(2,nodes(e,4))=inforce(2,nodes(e,4))+ninforce(8,1);
            end
        end
        nstress=(stab\glstress)';%得到物理坐标下节点的应力
        for i=1:4
        stress(:,nodes(e,i))=stress(:,nodes(e,i))+nstress(:,i);%节点应力求和
        end
    end
    stress=stress./numnd;%节点应力求均值
    delt=mindelt;
    a=exforce-inforce-c.*vel;
    a=a./mass;%得到节点加速度
    
    for n=1:numnode %添加加速度边界条件
        if(rem(n,(numelex+1))==1)
            a(1,n)=0;
            a(2,n)=0;
        end
    end 

    vel=vel+a*delt; %得到节点的速度和位移
    disp=disp+vel*delt;
    time(step)=step*delt;
    nxdisp(1,step)=disp(1,1326);
    nydisp(1,step)=disp(2,1326);
    nxdisp(2,step)=disp(1,663);
    nydisp(2,step)=disp(2,663);
    nxdisp(3,step)=disp(1,51);
    nydisp(3,step)=disp(2,51);
 if mod(step,10)==0
    numname=num2str(fix(step/10));
    name=['result',numname];
    name=strcat(name,'.plt');
    cd('result');
    fid_out=fopen(name,'w');
    fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
    fprintf(fid_out,'VARIABLES="x"  "y"  "u"  "v"  "sigmax"  "sigmay"  "sigmaxy"\n');
    fprintf(fid_out,'ZONE T="flow-field", N=%8d, E=%8d, ET=QUADRILATERAL, F=FEPOINT\n',numnode,numele);
    for i=1:numnode
       fprintf(fid_out,'%16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e\n',coornode(1,i),coornode(2,i),disp(1,i),disp(2,i),stress(1,i),stress(2,i),stress(3,i));
    end
    for i=1:numele
       fprintf(fid_out,'%8d %8d %8d %8d\n',nodes(i,1),nodes(i,2),nodes(i,3),nodes(i,4));
    end
    fclose(fid_out);
    cd('..');
   
 end
   
    step
    delt
end
plot(time,nxdisp(1,:),time,nydisp(1,:));
plot(time,nxdisp(2,:),time,nydisp(2,:));
plot(time,nxdisp(3,:),time,nydisp(3,:));

    










=======

%本程序采用的是CPS4单元即四节点四边形单元，计算格式采用的是显式格式，计算方法采用全积分，
%即在每个时间步都利用控制方程计算出所有节点的加速度，进而得到每个时间步节点的速度、位移和应力。
clc;
%前处理程序，直接建立模型，得到节点信息，单元信息等
lengthx=4;%模型大小
lengthy=2;
emodule=5e+6;%弹性模量
possion=0.25;%泊松比
density=0.0005;%密度
c=0.0;%阻尼
numelex=50;%x、y方向单元数
numeley=25;
numele=numelex*numeley;
numnode=(numelex+1)*(numeley+1);
nglx=2;%x、y方向高斯积分点数
ngly=2;
nnel=4;
coornode=[zeros(2,numnode)];%节点坐标
for i=1:numelex+1
    for j=1:numeley+1
        coornode(1,(j-1)*(numelex+1)+i)=(i-1)*lengthx/numelex;
        coornode(2,(j-1)*(numelex+1)+i)=(j-1)*lengthy/numeley;
    end
end

%单元对应节点坐标编号
nodes=[zeros(numele,4)];
mat1=1:numnode;
mat2=reshape(mat1,numelex+1,numeley+1);
for i=1:numele
    h=ceil(i/numelex);
    l=rem(i,numelex);
    if(l==0)
        l=numelex;
    end
    mat3=mat2(l:l+1,h:h+1);
    mat4=mat3(:);
    nodes(i,1:4)=[mat4([4,3,1,2])'];
end

%初始化节点内外力，施加力边界条件
exforce=[zeros(2,numnode)];
 %在最右端所有节点施加y方向大小为-1的集中力 
 for i=1:numeley+1
     exforce(2,(numelex+1)*i)=-1;
 end
inforce=[zeros(2,numnode)];
stress=[zeros(3,numnode)];

%计算每个节点质量，利用的一致质量矩阵原理，即每个单元的质量均分给它的节点
mass1=[zeros(1,numnode)];
for e=1:numele
    x1=coornode(1,nodes(e,1));
    x2=coornode(1,nodes(e,2));
    x3=coornode(1,nodes(e,3));
    x4=coornode(1,nodes(e,4));
    y1=coornode(2,nodes(e,1));
    y2=coornode(2,nodes(e,2));
    y3=coornode(2,nodes(e,3));
    y4=coornode(2,nodes(e,4));
    A=abs(0.5*(x2*y3+x3*y1+x1*y2-x3*y2-x1*y3-x2*y1));%单元面积
    A=A+abs(0.5*(x3*y4+x4*y2+x2*y3-x4*y3-x2*y4-x3*y2));
    m=density*A/4;
    mass1(1,nodes(e,1))=mass1(1,nodes(e,1))+m;
    mass1(1,nodes(e,2))=mass1(1,nodes(e,2))+m;
    mass1(1,nodes(e,3))=mass1(1,nodes(e,3))+m;
    mass1(1,nodes(e,4))=mass1(1,nodes(e,4))+m;
end
mass=[mass1;mass1];%节点质量矩阵

%初始化节点的位移、速度、加速度。
disp=[zeros(2,numnode)];
vel=[zeros(2,numnode)];
a=[zeros(2,numnode)];

%设置时间步
numstep=1000;
time=[zeros(1,numstep)];
D=(emodule/(1-possion^2))*[1,possion,0;possion,1,0;0,0,(1-possion)/2];%材料矩阵（平面应力类型）
[point2,weight2]=glqd2(nglx,ngly);%得到高斯积分点位置矩阵和权系数矩阵
stab=zeros(4,4);
nnd=zeros(1,numnode);%每个节点被单元共享次数
  for i=1:numele
      for j=1:4
          nnd(nodes(i,j))=nnd(nodes(i,j))+1;
      end
  end
  numnd=[nnd;nnd;nnd];
  nxdisp=zeros(3,numstep);
  nydisp=zeros(3,numstep);
 
  %开始时间步循环
for step=1:numstep
    inforce=[zeros(2,numnode)];
    stress=[zeros(3,numnode)];
    mindelt=0.0001;%初始设置单个时间步持续时间
    for e=1:numele %每个时间步单元循环
        glstress=[zeros(4,3)];
        X=zeros(1,4);
        Y=zeros(1,4);
        for i=1:4
         X(i)=coornode(1,nodes(e,i))+disp(1,nodes(e,i));
         Y(i)=coornode(2,nodes(e,i))+disp(2,nodes(e,i));
        end
         [a1,b1,c1,d1,angleA,angleB,angleC,angleD]=angle(X,Y);%得到每个单元内角
         %判断单元是否失效
       if angleA >=pi || angleB >=pi ||angleC >=pi || angleD >=pi 
          temp='quadrilateral element is concave - FEM can not solve';
          break;
       end
       
       %根据每个单元的最小边长计算最小时间步持续时间
       delt=min([a1,b1,c1,d1]);
       delt=0.7*delt/sqrt(emodule/density);
        if (delt<mindelt)
            mindelt=delt;
        end
        
        pk=0;
        for  inty=1:ngly %积分点循环
                y=point2(inty,2);      
                wty=weight2(inty,2);     
            for intx=1:nglx
             x=point2(intx,1);         
            wtx=weight2(intx,1);  
                [shape,dhdr,dhds]=glsb(x,y);%得到自然坐标下每个积分点的形函数及其偏导数
                jacob2=jacob(nnel,dhdr,dhds,X,Y); %得到高斯点的雅各比矩阵
                detjacob=det(jacob2);                        
                invjacob=inv(jacob2);
                [dhdx,dhdy]=glb(nnel,dhdr,dhds,invjacob);%得到物理坐标下节点形函数及其偏导数
                glB=feglB(nnel,dhdx,dhdy);%组装物理坐标下节点的B矩阵（应变矩阵）
                q=[zeros(8,1)];
                q(1,1)=disp(1,nodes(e,1));
                q(3,1)=disp(1,nodes(e,2));
                q(5,1)=disp(1,nodes(e,3));
                q(7,1)=disp(1,nodes(e,4));
                q(2,1)=disp(2,nodes(e,1));
                q(4,1)=disp(2,nodes(e,2));
                q(6,1)=disp(2,nodes(e,3));
                q(8,1)=disp(2,nodes(e,4));
                glnstrain=glB*q;%得到自然坐标的节点应变
                glnstress=D*glnstrain;%得到自然坐标的节点应力
                 pk=pk+1;
                 stab(pk,:)=shape;%组装应力转换矩阵
                glstress(pk,:)=glnstress';
                ninforce=detjacob*glB'*glnstress*wtx*wty;%得到物理坐标下节点的内力
                inforce(1,nodes(e,1))=inforce(1,nodes(e,1))+ninforce(1,1);
                inforce(2,nodes(e,1))=inforce(2,nodes(e,1))+ninforce(2,1);
                inforce(1,nodes(e,2))=inforce(1,nodes(e,2))+ninforce(3,1);
                inforce(2,nodes(e,2))=inforce(2,nodes(e,2))+ninforce(4,1);
                inforce(1,nodes(e,3))=inforce(1,nodes(e,3))+ninforce(5,1);
                inforce(2,nodes(e,3))=inforce(2,nodes(e,3))+ninforce(6,1);
                inforce(1,nodes(e,4))=inforce(1,nodes(e,4))+ninforce(7,1);
                inforce(2,nodes(e,4))=inforce(2,nodes(e,4))+ninforce(8,1);
            end
        end
        nstress=(stab\glstress)';%得到物理坐标下节点的应力
        for i=1:4
        stress(:,nodes(e,i))=stress(:,nodes(e,i))+nstress(:,i);%节点应力求和
        end
    end
    stress=stress./numnd;%节点应力求均值
    delt=mindelt;
    a=exforce-inforce-c.*vel;
    a=a./mass;%得到节点加速度
    
    for n=1:numnode %添加加速度边界条件
        if(rem(n,(numelex+1))==1)
            a(1,n)=0;
            a(2,n)=0;
        end
    end 

    vel=vel+a*delt; %得到节点的速度和位移
    disp=disp+vel*delt;
    time(step)=step*delt;
    nxdisp(1,step)=disp(1,1326);
    nydisp(1,step)=disp(2,1326);
    nxdisp(2,step)=disp(1,663);
    nydisp(2,step)=disp(2,663);
    nxdisp(3,step)=disp(1,51);
    nydisp(3,step)=disp(2,51);
 if mod(step,10)==0
    numname=num2str(fix(step/10));
    name=['result',numname];
    name=strcat(name,'.plt');
    cd('result');
    fid_out=fopen(name,'w');
    fprintf(fid_out,'TITLE="test case governed by poisson equation"\n');
    fprintf(fid_out,'VARIABLES="x"  "y"  "u"  "v"  "sigmax"  "sigmay"  "sigmaxy"\n');
    fprintf(fid_out,'ZONE T="flow-field", N=%8d, E=%8d, ET=QUADRILATERAL, F=FEPOINT\n',numnode,numele);
    for i=1:numnode
       fprintf(fid_out,'%16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e\n',coornode(1,i),coornode(2,i),disp(1,i),disp(2,i),stress(1,i),stress(2,i),stress(3,i));
    end
    for i=1:numele
       fprintf(fid_out,'%8d %8d %8d %8d\n',nodes(i,1),nodes(i,2),nodes(i,3),nodes(i,4));
    end
    fclose(fid_out);
    cd('..');
   
 end
   
    step
    delt
end
plot(time,nxdisp(1,:),time,nydisp(1,:));
plot(time,nxdisp(2,:),time,nydisp(2,:));
plot(time,nxdisp(3,:),time,nydisp(3,:));

    










>>>>>>> 0c2441c35348cc2cb0d4bfa4f93dc736c7d002fd
