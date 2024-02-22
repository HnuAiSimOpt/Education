%����Ԫ����ҵ���������嵥Ԫ�������
%ѧ��:S230200217 ����:����
clear
close all
clc
%---------------------
%  ����ڵ�͵�Ԫ��Ϣ
%---------------------
node=dlmread('node.txt',',');%����ڵ���Ϣ
elements=dlmread('elements.txt',',');%���뵥Ԫ��Ϣ
%---------------
%  ������Ʋ���
%---------------
x0=node(:,2:4);%�ڵ�����
nodes=elements(:,2:9);%��Ԫ����
nnode=size(x0,1);%�ڵ�����
nel=size(nodes,1);%��Ԫ����
nnel=size(nodes,2);%ÿ����Ԫ�ڵ���
ndof=3;%ÿ���ڵ����ɶ���
sdof=nnode*ndof;%ϵͳ�����ɶ���
edof=nnel*ndof;%ÿ����Ԫ���ɶ���
elastic=2.1e11;%����ģ��
poisson=0.3;%���ɱ�
ng=2;%��˹���ֵ���
%-------------------
%  ��ʼ�����������
%-------------------
F=zeros(sdof,1);%��ʼ���ڵ���
kk=sparse(sdof,sdof);%����նȾ���Ŀվ���
stress=zeros(nnode,7);%��ʼ��Ӧ������
ninode=zeros(nnode,1);%��ʼ���ڵ��������
mistress=zeros(nnel,1);%��ʼ��MisesӦ������
X=zeros(nnel,1);%��ʼ��X�������
Y=zeros(nnel,1);%��ʼ��Y�������
Z=zeros(nnel,1);%��ʼ��Z�������
%-----------------
%  ������������
%-----------------
D=Dc(4,elastic,poisson);%���Ծ���
[p,w]=gaussquadrature(ng);%��˹���ֵ�����ֵ�ͼ�Ȩϵ��
F(3*1160-1)=1000;%1160����y������
bcdof=[];%����Լ������Ŀվ���
for ibcdof=1:nnode
    if x0(ibcdof,3)==0
        bcdof=[bcdof ibcdof];%Լ������
    end
end
nbc=size(bcdof,2);%Լ������
for ini=1:nel
    for inj=1:nnel
        ninode(nodes(ini,inj))=ninode(nodes(ini,inj))+1;%ÿ���ڵ��õĴ���
    end
end
%---------------------
%  ����նȾ�����װ
%---------------------
for iel=1:nel
    nd=nodes(iel,:);%��iel����Ԫ�����ĸ��ڵ�
    for i=1:nnel
        X(i)=x0(nd(i),1);
        Y(i)=x0(nd(i),2);
        Z(i)=x0(nd(i),3);%��i����Ԫ�ĸ��ڵ�����
    end
    k=zeros(edof,edof);%��ʼ����Ԫ�նȾ���
    %--------------------------
    %  �ø�˹������Ԫ�նȾ���
    %--------------------------
    for iz=1:ng
        z=p(iz);
        wz=w(iz);
        for iy=1:ng
            y=p(iy);
            wy=w(iy);
            for ix=1:ng
                x=p(ix);
                wx=w(ix);
                [N,dndr,dnds,dndt]=shape(x,y,z);%�κ���
                J=Jacobain(nnel,dndr,dnds,dndt,X,Y,Z);%�ſ˱Ⱦ���
                B=Bc(J,nnel,dndr,dnds,dndt);%��ԪӦ�����
                k=k+(B')*D*B*wx*wy*wz*det(J);%���㵥Ԫ�նȾ���
            end
        end
    end
    %---------------
    %  ��װ�նȾ���
    %---------------
    kk=zk(kk,nd,k,nnel,ndof);%����նȾ���
end
%---------------
%  ʩ�ӱ߽�����
%---------------
kk=boundary(bcdof,nbc,kk,ndof);
%---------------
%  ����ڵ�λ��
%---------------
disp=kk\F;
%---------------------
%  ����Ӧ��������װ
%---------------------
for iel=1:nel
    nd=nodes(iel,:);%��iel����Ԫ�����ĸ��ڵ�
    for i=1:nnel
        X(i)=x0(nd(i),1);
        Y(i)=x0(nd(i),2);
        Z(i)=x0(nd(i),3);%��i����Ԫ�ĸ��ڵ�����
    end
    gstress=[];%��ʼ����˹��Ӧ������
    %----------------------------------
    %  �ø�˹���ֵ����˹���ֵ�Ӧ������
    %----------------------------------
    for iz=1:ng
        z=p(iz);
        for iy=1:ng
            y=p(iy);
            for ix=1:ng
                x=p(ix);
                [N,dndr,dnds,dndt]=shape(x,y,z);%�κ���
                J=Jacobain(nnel,dndr,dnds,dndt,X,Y,Z);%�ſ˱Ⱦ���
                B=Bc(J,nnel,dndr,dnds,dndt);%��ԪӦ�����
                stress_g=ssc(disp,ndof,nd,B,D);%��˹��Ӧ��
                gstress=[gstress stress_g];%��˹��Ӧ������
            end
        end
    end
    %---------------------
    %  �����˹���ֵ�Ӧ��
    %---------------------
    gstress=[gstress(:,5),gstress(:,6),gstress(:,8),gstress(:,7),...
        gstress(:,1),gstress(:,2),gstress(:,4),gstress(:,3)];
    %-----------
    %  Ӧ���ָ�
    %-----------
    a=sqrt(3);
    h=[recovery(-a,-a, a)
       recovery( a,-a, a)
       recovery( a, a, a)
       recovery(-a, a, a)
       recovery(-a,-a,-a)
       recovery( a,-a,-a)
       recovery( a, a,-a)
       recovery(-a, a,-a)];%Ӧ���ָ�����
    estress=h*transpose(gstress);%��Ԫ�ڵ�Ӧ��
    %---------------------------------
    %  ����MisesӦ�������鵥Ԫ�ڵ�Ӧ��
    %---------------------------------
    for mi=1:nnel
        s1=estress(mi,1);
        s2=estress(mi,2);
        s3=estress(mi,3);
        mistress(mi)=(((s1-s2)^2+(s2-s3)^2+(s3-s1)^2)/2)^0.5;%MisesӦ��
    end
    estress=[estress,mistress];%���鵥Ԫ�ڵ�Ӧ��
    %-----------------
    %  ��װ�ڵ���Ӧ��
    %-----------------
    for is=1:nnel
        stress(nd(is),:)=stress(nd(is),:)+estress(is,:);%���е�Ԫ���ӵĽڵ�Ӧ��
    end
end
%------------
%  Ӧ��ƽ��
%------------
for inode=1:nnode
    an=ninode(inode);
    stress(inode,:)=stress(inode,:)/an;%ÿ���ڵ�Ӧ��ƽ��
end
%--------------------------
%  ���������Ϊ�ļ�jieguo.plt��
%--------------------------
nel = size(nodes,1);
node= size(x0,1);
nodedata=fopen('jieguo.plt','w');  
fprintf(nodedata,'TITLE="data"\n');
fprintf(nodedata,'VARIABLES=,"X", "Y","Z" ,"sig1","sig2","sig3","sig4","sig5","sig6","sig7"\n');
fprintf(nodedata,'ZONE T="%d  "  ,  N=%d, E=%d, ET=BRICK, F=FEPOINT\n',1,node,nel);
sigm1=stress(:,1);
sigm2=stress(:,2);
sigm3=stress(:,3);
sigm4=stress(:,4);
sigm5=stress(:,5);
sigm6=stress(:,6);
sigm7=stress(:,7);
for i=1:node
    fprintf(nodedata,'%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f\n',...
      x0(i,1),x0(i,2),x0(i,3),sigm1(i),sigm2(i),sigm3(i),sigm4(i),sigm5(i),sigm6(i),sigm7(i));
end  
for i=1:nel
    for j=1:size(nodes,2)
        fprintf(nodedata,'%d       ',nodes(i,j));  
    end
    fprintf(nodedata,'\n');
end