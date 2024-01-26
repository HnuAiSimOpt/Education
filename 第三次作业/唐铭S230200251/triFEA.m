function [stress,u]=triFEA(nelx,nely,t,E,nu,F)



nodexnum=nelx+1;
nodeynum=nely+1;
nodenum=nodexnum*nodeynum;
elemnum=nelx*nely*2;
elemnode=zeros(elemnum,3);  %%%%���ǵ�Ԫ��Ӧ����Ϣ

elemnode(1:nely,1)=1:nely;
elemnode(1:nely,2)=elemnode(1:nely,1)+1;
elemnode(1:nely,3)=elemnode(1:nely,1)+nodeynum;
elemnode(nely+1:2*nely,1)=elemnode(1:nely,1)+1;
elemnode(nely+1:2*nely,2)=elemnode(1:nely,1)+nodeynum;
elemnode(nely+1:2*nely,3)=elemnode(nely+1:2*nely,2)+1; %%%��һ�е�Ԫ��Ӧ����Ϣ

      
for i=2:nelx
    elemnode(2*nely*(i-1)+1:2*nely*i,:)=elemnode(2*nely*(i-2)+1:2*nely*(i-1),:)+nodeynum;
end


nodeloc=zeros(nodenum,2);
for i=1:nodexnum
    nodeloc(nodeynum*(i-1)+1:nodeynum*i,1)=i-1;
    a=nodeynum-1:-1:0;
    nodeloc(nodeynum*(i-1)+1:nodeynum*i,2)=a;
end

K=zeros(2*nodenum,2*nodenum);
u=zeros(2*nodenum,1);
f=zeros(2*nodenum,1);
%% �߽�����

f(2*nodenum-2*nodeynum+1:2:2*nodenum-1,1)=F;%%%ʩ�Ӹ���
fixednode=[1:nodeynum*2]; %%%�̶�Լ��
allnode=[1:nodenum*2];  
freenode=setdiff(allnode,fixednode);

%% ���

for i=1:elemnum
    elemloc=zeros(3,2);  
    elemloc(1:3,:)=nodeloc(elemnode(i,1:3),:);
    ke=elemstiffness(E,nu,t,elemloc);
    K=stiffness(K,ke,elemnode(i,:));
end

u(freenode)=K(freenode,freenode)\f(freenode);
u(fixednode)=0;

stress=zeros(elemnum,3);

figure
for i=1:elemnum
 elemdisp=zeros(6,1);
 elemdisp(1:2:5)=u(2*elemnode(i,:)-1);
 elemdisp(2:2:6)=u(2*elemnode(i,:));
 elemloc=zeros(3,2);  
 elemloc(1:3,:)=nodeloc(elemnode(i,1:3),:);
 stress(i,:)=stresssolu(E,nu,elemloc,elemdisp);
 patch(elemloc(:,1),elemloc(:,2),stress(i,1));%%��ʾӦ��
end


%% ��ʾ���
nodexloc=nodeloc(:,1)+u(1:2:2*nodenum-1);
nodeyloc=nodeloc(:,2)+u(2:2:2*nodenum);
nodexlocold=nodeloc(:,1);
nodeylocold=nodeloc(:,2);

display_2D(elemnode,nodexloc,nodeyloc,nodexlocold,nodeylocold)  %%% ��ʾ���κ�ģ��
end