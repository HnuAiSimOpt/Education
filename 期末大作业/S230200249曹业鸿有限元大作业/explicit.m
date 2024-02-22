<<<<<<< HEAD

%��������õ���CPS4��Ԫ���Ľڵ��ı��ε�Ԫ�������ʽ���õ�����ʽ��ʽ�����㷽������ȫ���֣�
%����ÿ��ʱ�䲽�����ÿ��Ʒ��̼�������нڵ�ļ��ٶȣ������õ�ÿ��ʱ�䲽�ڵ���ٶȡ�λ�ƺ�Ӧ����
clc;
%ǰ�������ֱ�ӽ���ģ�ͣ��õ��ڵ���Ϣ����Ԫ��Ϣ��
lengthx=4;%ģ�ʹ�С
lengthy=2;
emodule=5e+6;%����ģ��
possion=0.25;%���ɱ�
density=0.0005;%�ܶ�
c=0.0;%����
numelex=50;%x��y����Ԫ��
numeley=25;
numele=numelex*numeley;
numnode=(numelex+1)*(numeley+1);
nglx=2;%x��y�����˹���ֵ���
ngly=2;
nnel=4;
coornode=[zeros(2,numnode)];%�ڵ�����
for i=1:numelex+1
    for j=1:numeley+1
        coornode(1,(j-1)*(numelex+1)+i)=(i-1)*lengthx/numelex;
        coornode(2,(j-1)*(numelex+1)+i)=(j-1)*lengthy/numeley;
    end
end

%��Ԫ��Ӧ�ڵ�������
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

%��ʼ���ڵ���������ʩ�����߽�����
exforce=[zeros(2,numnode)];
 %�����Ҷ����нڵ�ʩ��y�����СΪ-1�ļ����� 
 for i=1:numeley+1
     exforce(2,(numelex+1)*i)=-1;
 end
inforce=[zeros(2,numnode)];
stress=[zeros(3,numnode)];

%����ÿ���ڵ����������õ�һ����������ԭ����ÿ����Ԫ���������ָ����Ľڵ�
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
    A=abs(0.5*(x2*y3+x3*y1+x1*y2-x3*y2-x1*y3-x2*y1));%��Ԫ���
    A=A+abs(0.5*(x3*y4+x4*y2+x2*y3-x4*y3-x2*y4-x3*y2));
    m=density*A/4;
    mass1(1,nodes(e,1))=mass1(1,nodes(e,1))+m;
    mass1(1,nodes(e,2))=mass1(1,nodes(e,2))+m;
    mass1(1,nodes(e,3))=mass1(1,nodes(e,3))+m;
    mass1(1,nodes(e,4))=mass1(1,nodes(e,4))+m;
end
mass=[mass1;mass1];%�ڵ���������

%��ʼ���ڵ��λ�ơ��ٶȡ����ٶȡ�
disp=[zeros(2,numnode)];
vel=[zeros(2,numnode)];
a=[zeros(2,numnode)];

%����ʱ�䲽
numstep=1000;
time=[zeros(1,numstep)];
D=(emodule/(1-possion^2))*[1,possion,0;possion,1,0;0,0,(1-possion)/2];%���Ͼ���ƽ��Ӧ�����ͣ�
[point2,weight2]=glqd2(nglx,ngly);%�õ���˹���ֵ�λ�þ����Ȩϵ������
stab=zeros(4,4);
nnd=zeros(1,numnode);%ÿ���ڵ㱻��Ԫ�������
  for i=1:numele
      for j=1:4
          nnd(nodes(i,j))=nnd(nodes(i,j))+1;
      end
  end
  numnd=[nnd;nnd;nnd];
  nxdisp=zeros(3,numstep);
  nydisp=zeros(3,numstep);
 
  %��ʼʱ�䲽ѭ��
for step=1:numstep
    inforce=[zeros(2,numnode)];
    stress=[zeros(3,numnode)];
    mindelt=0.0001;%��ʼ���õ���ʱ�䲽����ʱ��
    for e=1:numele %ÿ��ʱ�䲽��Ԫѭ��
        glstress=[zeros(4,3)];
        X=zeros(1,4);
        Y=zeros(1,4);
        for i=1:4
         X(i)=coornode(1,nodes(e,i))+disp(1,nodes(e,i));
         Y(i)=coornode(2,nodes(e,i))+disp(2,nodes(e,i));
        end
         [a1,b1,c1,d1,angleA,angleB,angleC,angleD]=angle(X,Y);%�õ�ÿ����Ԫ�ڽ�
         %�жϵ�Ԫ�Ƿ�ʧЧ
       if angleA >=pi || angleB >=pi ||angleC >=pi || angleD >=pi 
          temp='quadrilateral element is concave - FEM can not solve';
          break;
       end
       
       %����ÿ����Ԫ����С�߳�������Сʱ�䲽����ʱ��
       delt=min([a1,b1,c1,d1]);
       delt=0.7*delt/sqrt(emodule/density);
        if (delt<mindelt)
            mindelt=delt;
        end
        
        pk=0;
        for  inty=1:ngly %���ֵ�ѭ��
                y=point2(inty,2);      
                wty=weight2(inty,2);     
            for intx=1:nglx
             x=point2(intx,1);         
            wtx=weight2(intx,1);  
                [shape,dhdr,dhds]=glsb(x,y);%�õ���Ȼ������ÿ�����ֵ���κ�������ƫ����
                jacob2=jacob(nnel,dhdr,dhds,X,Y); %�õ���˹����Ÿ��Ⱦ���
                detjacob=det(jacob2);                        
                invjacob=inv(jacob2);
                [dhdx,dhdy]=glb(nnel,dhdr,dhds,invjacob);%�õ����������½ڵ��κ�������ƫ����
                glB=feglB(nnel,dhdx,dhdy);%��װ���������½ڵ��B����Ӧ�����
                q=[zeros(8,1)];
                q(1,1)=disp(1,nodes(e,1));
                q(3,1)=disp(1,nodes(e,2));
                q(5,1)=disp(1,nodes(e,3));
                q(7,1)=disp(1,nodes(e,4));
                q(2,1)=disp(2,nodes(e,1));
                q(4,1)=disp(2,nodes(e,2));
                q(6,1)=disp(2,nodes(e,3));
                q(8,1)=disp(2,nodes(e,4));
                glnstrain=glB*q;%�õ���Ȼ����Ľڵ�Ӧ��
                glnstress=D*glnstrain;%�õ���Ȼ����Ľڵ�Ӧ��
                 pk=pk+1;
                 stab(pk,:)=shape;%��װӦ��ת������
                glstress(pk,:)=glnstress';
                ninforce=detjacob*glB'*glnstress*wtx*wty;%�õ����������½ڵ������
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
        nstress=(stab\glstress)';%�õ����������½ڵ��Ӧ��
        for i=1:4
        stress(:,nodes(e,i))=stress(:,nodes(e,i))+nstress(:,i);%�ڵ�Ӧ�����
        end
    end
    stress=stress./numnd;%�ڵ�Ӧ�����ֵ
    delt=mindelt;
    a=exforce-inforce-c.*vel;
    a=a./mass;%�õ��ڵ���ٶ�
    
    for n=1:numnode %��Ӽ��ٶȱ߽�����
        if(rem(n,(numelex+1))==1)
            a(1,n)=0;
            a(2,n)=0;
        end
    end 

    vel=vel+a*delt; %�õ��ڵ���ٶȺ�λ��
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

%��������õ���CPS4��Ԫ���Ľڵ��ı��ε�Ԫ�������ʽ���õ�����ʽ��ʽ�����㷽������ȫ���֣�
%����ÿ��ʱ�䲽�����ÿ��Ʒ��̼�������нڵ�ļ��ٶȣ������õ�ÿ��ʱ�䲽�ڵ���ٶȡ�λ�ƺ�Ӧ����
clc;
%ǰ�������ֱ�ӽ���ģ�ͣ��õ��ڵ���Ϣ����Ԫ��Ϣ��
lengthx=4;%ģ�ʹ�С
lengthy=2;
emodule=5e+6;%����ģ��
possion=0.25;%���ɱ�
density=0.0005;%�ܶ�
c=0.0;%����
numelex=50;%x��y����Ԫ��
numeley=25;
numele=numelex*numeley;
numnode=(numelex+1)*(numeley+1);
nglx=2;%x��y�����˹���ֵ���
ngly=2;
nnel=4;
coornode=[zeros(2,numnode)];%�ڵ�����
for i=1:numelex+1
    for j=1:numeley+1
        coornode(1,(j-1)*(numelex+1)+i)=(i-1)*lengthx/numelex;
        coornode(2,(j-1)*(numelex+1)+i)=(j-1)*lengthy/numeley;
    end
end

%��Ԫ��Ӧ�ڵ�������
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

%��ʼ���ڵ���������ʩ�����߽�����
exforce=[zeros(2,numnode)];
 %�����Ҷ����нڵ�ʩ��y�����СΪ-1�ļ����� 
 for i=1:numeley+1
     exforce(2,(numelex+1)*i)=-1;
 end
inforce=[zeros(2,numnode)];
stress=[zeros(3,numnode)];

%����ÿ���ڵ����������õ�һ����������ԭ����ÿ����Ԫ���������ָ����Ľڵ�
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
    A=abs(0.5*(x2*y3+x3*y1+x1*y2-x3*y2-x1*y3-x2*y1));%��Ԫ���
    A=A+abs(0.5*(x3*y4+x4*y2+x2*y3-x4*y3-x2*y4-x3*y2));
    m=density*A/4;
    mass1(1,nodes(e,1))=mass1(1,nodes(e,1))+m;
    mass1(1,nodes(e,2))=mass1(1,nodes(e,2))+m;
    mass1(1,nodes(e,3))=mass1(1,nodes(e,3))+m;
    mass1(1,nodes(e,4))=mass1(1,nodes(e,4))+m;
end
mass=[mass1;mass1];%�ڵ���������

%��ʼ���ڵ��λ�ơ��ٶȡ����ٶȡ�
disp=[zeros(2,numnode)];
vel=[zeros(2,numnode)];
a=[zeros(2,numnode)];

%����ʱ�䲽
numstep=1000;
time=[zeros(1,numstep)];
D=(emodule/(1-possion^2))*[1,possion,0;possion,1,0;0,0,(1-possion)/2];%���Ͼ���ƽ��Ӧ�����ͣ�
[point2,weight2]=glqd2(nglx,ngly);%�õ���˹���ֵ�λ�þ����Ȩϵ������
stab=zeros(4,4);
nnd=zeros(1,numnode);%ÿ���ڵ㱻��Ԫ�������
  for i=1:numele
      for j=1:4
          nnd(nodes(i,j))=nnd(nodes(i,j))+1;
      end
  end
  numnd=[nnd;nnd;nnd];
  nxdisp=zeros(3,numstep);
  nydisp=zeros(3,numstep);
 
  %��ʼʱ�䲽ѭ��
for step=1:numstep
    inforce=[zeros(2,numnode)];
    stress=[zeros(3,numnode)];
    mindelt=0.0001;%��ʼ���õ���ʱ�䲽����ʱ��
    for e=1:numele %ÿ��ʱ�䲽��Ԫѭ��
        glstress=[zeros(4,3)];
        X=zeros(1,4);
        Y=zeros(1,4);
        for i=1:4
         X(i)=coornode(1,nodes(e,i))+disp(1,nodes(e,i));
         Y(i)=coornode(2,nodes(e,i))+disp(2,nodes(e,i));
        end
         [a1,b1,c1,d1,angleA,angleB,angleC,angleD]=angle(X,Y);%�õ�ÿ����Ԫ�ڽ�
         %�жϵ�Ԫ�Ƿ�ʧЧ
       if angleA >=pi || angleB >=pi ||angleC >=pi || angleD >=pi 
          temp='quadrilateral element is concave - FEM can not solve';
          break;
       end
       
       %����ÿ����Ԫ����С�߳�������Сʱ�䲽����ʱ��
       delt=min([a1,b1,c1,d1]);
       delt=0.7*delt/sqrt(emodule/density);
        if (delt<mindelt)
            mindelt=delt;
        end
        
        pk=0;
        for  inty=1:ngly %���ֵ�ѭ��
                y=point2(inty,2);      
                wty=weight2(inty,2);     
            for intx=1:nglx
             x=point2(intx,1);         
            wtx=weight2(intx,1);  
                [shape,dhdr,dhds]=glsb(x,y);%�õ���Ȼ������ÿ�����ֵ���κ�������ƫ����
                jacob2=jacob(nnel,dhdr,dhds,X,Y); %�õ���˹����Ÿ��Ⱦ���
                detjacob=det(jacob2);                        
                invjacob=inv(jacob2);
                [dhdx,dhdy]=glb(nnel,dhdr,dhds,invjacob);%�õ����������½ڵ��κ�������ƫ����
                glB=feglB(nnel,dhdx,dhdy);%��װ���������½ڵ��B����Ӧ�����
                q=[zeros(8,1)];
                q(1,1)=disp(1,nodes(e,1));
                q(3,1)=disp(1,nodes(e,2));
                q(5,1)=disp(1,nodes(e,3));
                q(7,1)=disp(1,nodes(e,4));
                q(2,1)=disp(2,nodes(e,1));
                q(4,1)=disp(2,nodes(e,2));
                q(6,1)=disp(2,nodes(e,3));
                q(8,1)=disp(2,nodes(e,4));
                glnstrain=glB*q;%�õ���Ȼ����Ľڵ�Ӧ��
                glnstress=D*glnstrain;%�õ���Ȼ����Ľڵ�Ӧ��
                 pk=pk+1;
                 stab(pk,:)=shape;%��װӦ��ת������
                glstress(pk,:)=glnstress';
                ninforce=detjacob*glB'*glnstress*wtx*wty;%�õ����������½ڵ������
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
        nstress=(stab\glstress)';%�õ����������½ڵ��Ӧ��
        for i=1:4
        stress(:,nodes(e,i))=stress(:,nodes(e,i))+nstress(:,i);%�ڵ�Ӧ�����
        end
    end
    stress=stress./numnd;%�ڵ�Ӧ�����ֵ
    delt=mindelt;
    a=exforce-inforce-c.*vel;
    a=a./mass;%�õ��ڵ���ٶ�
    
    for n=1:numnode %��Ӽ��ٶȱ߽�����
        if(rem(n,(numelex+1))==1)
            a(1,n)=0;
            a(2,n)=0;
        end
    end 

    vel=vel+a*delt; %�õ��ڵ���ٶȺ�λ��
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
