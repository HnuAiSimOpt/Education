%% ע��
%{                 
                                                    ����Ԫ��ĩ����ҵ
  �������ξ���      ѧ�ţ�S230200192
  1.����Ԫ���򣺽��ƽ��Ӧ������
  2.��Ԫ���ͣ�4�ڵ��ı��ε�Ԫ
  3.�غɱ�ʩ��������һ�࣬���෴�����ǹ̶���
%}
clear
format short
first_time=cputime;             % ���ڼ������ִ����ʱ��
%% ��ʼ��
% �������
E=2.1e5;                        %����ģ��
miu=0.3;                        %���ɱ�
% ����ڵ���Ϣ����
node=importdata('nodes1.txt');  % �ڵ�������Ϣ����ʾ�ڵ�x��y���������
ele=importdata('element1.txt'); % ��Ԫ�ڵ���Ϣ����ʾÿ����Ԫ�ϵĽڵ����
nel=length(ele(:,1));           % ��Ԫ�ܸ���
nnode=length(node(:,1));        % �ڵ��ܸ���
nnel=4;                         % ÿ��Ԫ�ڵ���
ndof=2;                         % ÿ�ڵ����ɶ���
sdof=nnode*ndof;                % �����ɶ���
edof=nnel*ndof;                 % ÿ��Ԫ���ɶ���
% ��ʾ���񻮷����
figure(1)
hold on
axis off
axis equal
xp=[];
yp=[];
for ie=1:nel
    for j=1:nnel+1
        j1=mod(j-1,nnel)+1;
        xp(j)=node(ele(ie,j1),1);
        yp(j)=node(ele(ie,j1),2);

    end
    plot(xp,yp,'b-')
end
%% ʩ�ӱ߽��������غ�
fix_node=[26,  28, 12, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141,142, 143, 144, 145, 146];
bcdof=[];               % ��Լ�������ɶ�
bcval=[];               % ��Լ�����ɶȶ�Ӧ��λ��Լ��ֵ
for i=1:length(fix_node)             
    bcdof=[bcdof 2*fix_node(i)-1 2*fix_node(i)];
    bcval=[bcval 0 0];
end
% ʩ���غ�
ff(54,1)=-1000;                             % �����Ͻǽڵ�ʩ��һ��y��������µļ���������СΪ1000N

%% ��Ԫ�նȾ������������նȾ������װ
% �����ʼ��
ff=zeros(sdof,1);			% ϵͳ�غɳ�ʼ��
k=zeros(edof,edof);		    % ��Ԫ�նȳ�ʼ��
kk=sparse(sdof,sdof);		% ����նȳ�ʼ��
disp=zeros(sdof,1);         % ����λ�Ƴ�ʼ��
index=zeros(edof,1);		% ����Ԫ��8�����ɶȶ�Ӧ����������ɶȱ��
B=sparse(3,edof);           % Ӧ�����B
D=E/(1-miu*miu)*[1  miu 0;
   miu  1  0;
   0  0  (1-miu)/2];        % ƽ��Ӧ������ı�������D

G=2;                        % ��˹���ֻ��ֵ�ĸ�����G=1ʱΪһ���˹����
[point,weight]=fun_GaussQuadrature(G);
                            % fun_GaussQuadrature()�������ظ�˹���ֵĻ��ֵ������Ȩֵ�����庬�����Ӻ������н���
for iel=1:nel               % ��������Ԫ
    for i=1:nnel            % ����ÿ��Ԫ���ڵ�
        nd(i)=ele(iel,i);             % ���� nd ������ȡÿ����Ԫ�Ľڵ����
        xcoord(i)=node(nd(i),1);      % ���� xcoord ������ȡÿ����Ԫ�Ľڵ����������ռ��µ�x����
        ycoord(i)=node(nd(i),2);      % ���� ycoord ������ȡÿ����Ԫ�Ľڵ����������ռ��µ�y����
    end
    
    % ��ʼ����նȾ���
    for ix=1:2                      % ����x����ĸ����ֵ�
        x=point(ix,1);              % ��ȡx������ֵ������ֵ
        wx=weight(ix,1);            % ��ȡx������ֵ��Ȩ��
        for iy=1:2                  % ����y����ĸ����ֵ�
            y=point(iy,2);          % ��ȡy������ֵ������ֵ
            wy=weight(iy,2);        % ��ȡy������ֵ��Ȩ��
            
            [shape,dNdr,dNds]=fun_shapeFunction(x,y);             % ������ֵ���κ�����ֵ�Ͷ���Ȼ����r��s�ĵ���
            
            Jacobian=fun_Jacobian(nnel,dNdr,dNds,xcoord,ycoord);  % �����ſɱȾ���
            
            detJacobian=det(Jacobian);                            % �����ſɱȾ��������ʽ           
            
            invJacobian=inv(Jacobian);                            % �����ſɱȾ������
            
            [dNdx,dNdy]=fun_dNdx_dNdy(nnel,dNdr,dNds,invJacobian);% �����κ���������ռ�����x��y�ĵ���
                                                                  % dNdx��dNdy��4ά���������ֱ��ʾ��i���κ�����x��y�ĵ���
            
            B=fun_B(nnel,dNdx,dNdy);                              % ��װB����
   
            k=k+B'*D*B*wx*wy*detJacobian;                         % ��װ��Ԫ�նȾ���
            
        end
    end
    
    index=fun_index(nd,nnel,ndof);      % ��ȡ��Ԫ��ϵͳ���ɶȱ��
    
   	kk=fun_assembleKK(kk,k,index);      % ��װ����նȾ���
end

%% ��λ�Ʊ߽����������Խ���Ԫ�ظ�1��
[kk,ff]=fun_boundary(kk,ff,bcdof,bcval);

%% �������
[L,U]=lu(kk);
temp=L\ff;
disp=U\temp ;                               % �õ��ڵ�λ��

solve_time = cputime-first_time             % ������ʱ��

% ���κ������
disp=load('disp.txt');
disp_x=disp(1:2:sdof);
disp_y=disp(2:2:sdof);
dis=[disp_x,disp_y];
node_new=node+dis;
%% ��ԪӦ��Ӧ��Ӧ�����
const=[1.866 -0.5 0.134 -0.5;         % ������������Ӧ�����ƣ��ø�˹���Ӧ��ֵ����ֵ�ڵ��Ӧ��ֵ
      -0.5 1.866 -0.5 0.134;
      0.134 -0.5 1.866 -0.5;
      -0.5 0.134 -0.5 1.866];
 stress=zeros(nel,4,3);
 strain=zeros(nel,4,3);
for iel=1:nel                         % ��������Ԫ
    for i=1:nnel                      % ����ÿ��Ԫ���ڵ�
        nd(i)=ele(iel,i);             % ���� nd ������ȡÿ����Ԫ�Ľڵ����
        xcoord(i)=node(nd(i),1);      % ���� xcoord ������ȡÿ����Ԫ�Ľڵ����������ռ��µ�x����
        ycoord(i)=node(nd(i),2);      % ���� ycoord ������ȡÿ����Ԫ�Ľڵ����������ռ��µ�y����
    end
    pk=0;
    stress_Gauss=zeros(3,4);
    strain_Gauss=zeros(3,4);
    for ix=1:2                        % ����x����ĸ����ֵ�
        x=point(ix,1);                % ��ȡx������ֵ������ֵ
        wx=weight(ix,1);              % ��ȡx������ֵ��Ȩ��
        for iy=1:2                    % ����y����ĸ����ֵ�
            y=point(iy,2);            % ��ȡy������ֵ������ֵ
            wy=weight(iy,2);          % ��ȡy������ֵ��Ȩ��
            
            [shape,dNdr,dNds]=fun_shapeFunction(x,y);             % ������ֵ���κ�����ֵ�Ͷ���Ȼ����r��s�ĵ���
            
            Jacobian=fun_Jacobian(nnel,dNdr,dNds,xcoord,ycoord);  % �����ſɱȾ���
            
            detJacobian=det(Jacobian);                            % �����ſɱȾ��������ʽ           
            
            invJacobian=inv(Jacobian);                            % �����ſɱȾ������
            
            [dNdx,dNdy]=fun_dNdx_dNdy(nnel,dNdr,dNds,invJacobian);% �����κ���������ռ�����x��y�ĵ���
                                                                  % dNdx��dNdy��4ά���������ֱ��ʾ��i���κ�����x��y�ĵ���
            
            B=fun_B(nnel,dNdx,dNdy);                              % ��װB����
   
            index=fun_index(nd,nnel,ndof);                        % ��Ԫ�ڵ����ɶ�����
            for i=1:edof
                disp_ele(i,1)=disp(index(i));                     % ��Ԫ�ڵ�λ��
            end            
            strain_ele=B*disp_ele;                                % ���㵥ԪӦ��
            stress_ele=D*strain_ele;                              % ���㵥ԪӦ��
            pk=pk+1;
            stress_Gauss(:,pk)=stress_ele;
            strain_Gauss(:,pk)=strain_ele;
        end
    end
    for i=1:3
        for j=1:4
            for n=1:4
                stress(iel,j,i)= stress(iel,j,i)+const(j,n)*stress_Gauss(i,n);
                strain(iel,j,i)= strain(iel,j,i)+const(j,n)*strain_Gauss(i,n);
            end
        end
    end
    
end

% �ҳ�ÿ�ڵ����ڵĵ�Ԫ
neigh_node = cell(nnode,1);
neigh_node_ind = cell(nnode,1);
indneigh=zeros(1,nnode);
for i=1:nel
    for j=1:4
        indneigh(ele(i,j))=indneigh(ele(i,j))+1;
        neigh_node{ele(i,j)}(indneigh(ele(i,j)))=i;
        neigh_node_ind{ele(i,j)}(indneigh(ele(i,j)))=j;
    end
end

% ����ڵ㴦��Ӧ��Ӧ��
stress_node=zeros(3,nnode);
strain_node=zeros(3,nnode);
for inode=1:nnode
    numel= indneigh(inode);
    for i=1:numel
        ind_nel= neigh_node{inode}(i);
        ind_nod=neigh_node_ind{inode}(i);
        for j=1:3
            stress_node(j,inode)=stress_node(j,inode)+stress(ind_nel,ind_nod,j);
            strain_node(j,inode)=strain_node(j,inode)+strain(ind_nel,ind_nod,j);
        end
    end
    stress_node(:,inode)=stress_node(:,inode)/numel;       % �ڵ㴦��Ӧ��Ϊ��ͬ���ڵ�Ԫĥƽ��Ӧ����ƽ��ֵ
    strain_node(:,inode)=strain_node(:,inode)/numel;
end

%% ��ͼ
  fun_picture(node,disp,stress_node,strain_node);
%% ���
% ���������
fun_check(disp,stress_node,strain_node);    

% ������м�������report_all.txt�ļ�
%-------------------------------------------------------------
fid_out=fopen('report_all.txt','w');
fprintf(fid_out,'TITLE="���������"\n');
fprintf(fid_out,'VARIABLES="x" "y" "u" "v" "sigax"  "sigmay" "sigmaxy" "epsilonx"  "epsilony" "epsilonxy"\n');
fprintf(fid_out,' N= %8d,E=%8d,Element Type = QUADRILATERA\n',nnode,nel);
for i=1:nnode
    fprintf(fid_out,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',node(i,1),node(i,2),disp(2*i-1),disp(2*i),stress_node(1,i),stress_node(2,i),stress_node(3,i),strain_node(1,i),strain_node(2,i),strain_node(3,i));
end
for i=1:nel
    fprintf(fid_out,'%8d%8d%8d%8d\n',ele(i,1),ele(i,2),ele(i,3),ele(i,4));
end
%% �Ӻ���
% ��ȡ��Ԫ���������ɶȱ��
function index=fun_index(nd,nnel,ndof)
% nd����Ԫ�ڵ��
% nnel��ÿ��Ԫ�ڵ���
% ndof��ÿ�ڵ����ɶ�
 edof = nnel*ndof;
   k=0;
   for i=1:nnel
     start = (nd(i)-1)*ndof;
       for j=1:ndof
         k=k+1;
         index(k)=start+j;
       end
   end

end
% �����κ�������������x��y�ĵ���
function [dNdx,dNdy]=fun_dNdx_dNdy(nnel,dNdr,dNds,invJacobian) 
for i=1:nnel
 dNdx(i)=invJacobian(1,1)*dNdr(i)+invJacobian(1,2)*dNds(i);
 dNdy(i)=invJacobian(2,1)*dNdr(i)+invJacobian(2,2)*dNds(i);
end
end
% ���ظ�˹���ֵ㴦���κ���ֵ���κ�������Ȼ����ĵ���ֵ
function [shape,dNdr,dNds]=fun_shapeFunction(r,s)
% r��sΪ���ֵ�����ֵ
% ������Ȼ����ϵ��ÿ���ڵ���κ���
 shape(1)=0.25*(1-r)*(1-s);
 shape(2)=0.25*(1+r)*(1-s);
 shape(3)=0.25*(1+r)*(1+s);
 shape(4)=0.25*(1-r)*(1+s);
 
% �����κ�������Ȼ����r�ĵ���
 dNdr(1)=-0.25*(1-s);
 dNdr(2)=0.25*(1-s);
 dNdr(3)=0.25*(1+s);
 dNdr(4)=-0.25*(1+s);
 
% �����κ�������Ȼ����s�ĵ���
 dNds(1)=-0.25*(1-r);
 dNds(2)=-0.25*(1+r);
 dNds(3)=0.25*(1+r);
 dNds(4)=0.25*(1-r);
end
% �����ſɱȾ���
function Jacobian=fun_Jacobian(nnel,dNdr,dNds,xcoord,ycoord)
% nnel:��Ԫ�ڵ���
% dNdr��dNds���κ�������Ȼ����ĵ���
% xcoord,ycoord������ռ��µ�Ԫ�ĸ��ڵ������
 Jacobian=zeros(2,2);
 for i=1:nnel
     Jacobian(1,1)=Jacobian(1,1)+dNdr(i)*xcoord(i);
     Jacobian(1,2)=Jacobian(1,2)+dNdr(i)*ycoord(i);
     Jacobian(2,1)=Jacobian(2,1)+dNds(i)*xcoord(i);
     Jacobian(2,2)=Jacobian(2,2)+dNds(i)*ycoord(i);
 end
end
% �����ά��˹���ֵ�����ڵ�
function [point,weight]=fun_GaussQuadrature(i)
% point �� weight ����ĵ�һ�б�ʾx����ֵĽڵ��Ȩ�� 
% point �� weight ����ĵڶ��б�ʾy����ֵĽڵ��Ȩ��
% �м�������ڵ� point �� weight ���м���
   point=zeros(i,2);        % i��ʾ��i������ڵ㣬2ָ��x���y��
   weight=zeros(i,2);
   
   if i==1                  % 1���˹����
       point=[0,0];
       weight=[2,2];
   elseif i==2              % 2���˹����
       point=[-0.577,-0.577;
               0.577,0.577];
       weight=[1,1
               1,1];
   end
    
   
end
% ��װ����նȾ���
function kk=fun_assembleKK(kk,k,index)
edof = length(index);
for i=1:edof
    ii=index(i);
    for j=1:edof
        jj=index(j);
        kk(ii,jj)=kk(ii,jj)+k(i,j);
    end
end
end
% ��װӦ�����B
function B=fun_B(nnel,dNdx,dNdy) 
for i=1:nnel
    i1=(i-1)*2+1;  
    i2=i1+1;
    B(1,i1)=dNdx(i);
    B(2,i2)=dNdy(i);
    B(3,i1)=dNdy(i);
    B(3,i2)=dNdx(i);
end
end
% ����߽����������öԽ���Ԫ�ظ�1����ȥ�������Է������Լ����
function [kk,ff]=fun_boundary(kk,ff,bcdof,bcval)
 n=length(bcdof);
 sdof=length(kk);

 for i=1:n
    c=bcdof(i);
    for j=1:sdof
       kk(c,j)=0;
    end

    kk(c,c)=1;
    ff(c)=bcval(i);
 end
end
% ������ͼ
function fun_picture(node,disp,stress_node,strain_node)
x=node(:,1);
y=node(:,2);
% ����λ����ͼ
[m,n]=size(disp);
z=[];
for i=1:2:m-1
    z=[z;sqrt(disp(i)*disp(i)+disp(i+1)*disp(i+1))];
end
[X,Y,Z]=griddata(x,y,z,linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');  %��ֵ
%ȥ��ģ���ⲿ�ֵ�ֵ
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z(j,i)=NaN;
        end
    end
end
%���Ƶȸ���ͼ
figure('name','��λ����ͼ');
contourf(X,Y,Z,'LineStyle','none');
colormap jet
hold on;
colorbar;
axis equal
title('��λ����ͼ');

%����Ӧ����ͼ
[m1,n2]=size(stress_node);
z2=[];
for i=1:n2
        z2=[z2 sqrt(stress_node(1,i)*stress_node(1,i)+stress_node(2,i)*stress_node(2,i))];
end
[X,Y,Z0]=griddata(x,y,z2,linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear'); 
%ȥ��ģ���ⲿ�ֵ�ֵ
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z0(j,i)=NaN;
        end
    end
end
%���Ƶȸ���ͼ
figure('name','��Ӧ����ͼ');
contourf(X,Y,Z0,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('��Ӧ����ͼ');

%%x��Ӧ��
[X,Y,Z1]=griddata(x,y,stress_node(1,:),linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');
%ȥ��ģ���ⲿ�ֵ�ֵ
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z1(j,i)=NaN;
        end
    end
end
%���Ƶȸ���ͼ
figure('name','x�᷽��Ӧ����ͼ');
contourf(X,Y,Z1,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('x�᷽��Ӧ����ͼ');

%%y��Ӧ��
%��ֵ
[X,Y,Z2]=griddata(x,y,stress_node(2,:),linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');
%ȥ��ģ���ⲿ�ֵ�ֵ
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z2(j,i)=NaN;
        end
    end
end
%���Ƶȸ���ͼ
figure('name','y�᷽��Ӧ����ͼ');
contourf(X,Y,Z2,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('y�᷽��Ӧ����ͼ');

%%��Ӧ��
[X,Y,Z3]=griddata(x,y,stress_node(3,:),linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');
%ȥ��ģ���ⲿ�ֵ�ֵ
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z3(j,i)=NaN;
        end
    end
end
%���Ƶȸ���ͼ
figure('name','��Ӧ����ͼ');
contourf(X,Y,Z3,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('��Ӧ����ͼ');

%%����Ӧ����ͼ
%%x��Ӧ��
[X,Y,Z1]=griddata(x,y,strain_node(1,:),linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');
%ȥ��ģ���ⲿ�ֵ�ֵ
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z1(j,i)=NaN;
        end
    end
end
%���Ƶȸ���ͼ
figure('name','x�᷽��Ӧ����ͼ');
contourf(X,Y,Z1,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('x�᷽��Ӧ����ͼ');

%%y��Ӧ��
%��ֵ
[X,Y,Z2]=griddata(x,y,strain_node(2,:),linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');
%ȥ��ģ���ⲿ�ֵ�ֵ
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z2(j,i)=NaN;
        end
    end
end
%���Ƶȸ���ͼ
figure('name','y�᷽��Ӧ����ͼ');
contourf(X,Y,Z2,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('y�᷽��Ӧ����ͼ');

%%��Ӧ��
[X,Y,Z3]=griddata(x,y,strain_node(3,:),linspace(min(x),max(x),2000)',linspace(min(y),max(y),1000),'linear');
%ȥ��ģ���ⲿ�ֵ�ֵ
for i=1:2000
    for j=1:1000
        if sqrt((i-1000)^2+(j-500)^2)<200
            Z3(j,i)=NaN;
        end
    end
end
%���Ƶȸ���ͼ
figure('name','��Ӧ����ͼ');
contourf(X,Y,Z3,'LineStyle','none')  ;
colormap jet
hold on;
colorbar;
axis equal;
title('��Ӧ����ͼ');
end
% ������Ա�������
function fun_check(disp,stress,strain)
stress_node=stress';
strain_node=strain';
% ����λ�ơ�Ӧ����Ӧ�����
u0=load('u_abaqus.txt');      
[m,n]=size(u0);
for i=1:m
    u(2*i-1, 1)=u0(i, 1);
    u(2*i,1)=u0(i,2);
end
s= load('s_abaqus.txt');
E= load('E_abaqus.txt');
% ������� 
d_u=disp-u;         %λ�Ʋ�ֵ
d_s=stress_node-s;  %Ӧ����ֵ
d_E=strain_node-E;  %Ӧ���ֵ

for i=1:2*m
    eu(i,1)=abs(d_u(i,1)/disp(i,1));
end
for j=1:m
    esx(j,1)=abs(d_s(j,1)/stress_node(j,1));
    esy(j,1)=abs(d_s(j,2)/stress_node(j,2));
    esxy(j,1)=abs(d_s(j,3)/stress_node(j,3));
    
    eEx(j,1)=abs(d_E(j,1)/strain_node(j,1));
    eEy(j,1)=abs(d_E(j,2)/strain_node(j,2));
    eExy(j,1)=abs(d_E(j,3)/strain_node(j,3));
end
es=[esx;esy;esxy];
eE=[eEx;eEy;eExy];
% ��������ֵ
eu_total=mean(eu);
es_total=mean(es);
eE_total=mean(eE);
end