function [u]=squareFEM(nelx,nely,E,nu,F,t)
%%% 2D  4 nodes model
%% basic information
elemnum=nelx*nely;

xdot=nelx+1;
ydot=nely+1;
%nodenum=xdot*ydot;


nodenum=(2*xdot-1)*(2*ydot-1)-nelx*nely;

%% numbering nodes
num=zeros(elemnum,8);


    
    for i=1:nelx
        
        num((i-1)*nely+1:i*nely,1)=ydot*(i-1)+1:ydot*(i-1)+nely;
        num((i-1)*nely+1:i*nely,2)=num((i-1)*nely+1:i*nely,1)+1;
        num((i-1)*nely+1:i*nely,3)=num((i-1)*nely+1:i*nely,2)+ydot;
        num((i-1)*nely+1:i*nely,4)=num((i-1)*nely+1:i*nely,1)+ydot;
    
        %%%%%%
        num((i-1)*nely+1:i*nely,5)=nely*(i-1)+1+xdot*ydot:nely*(i-1)+nely+xdot*ydot;
        num((i-1)*nely+1:i*nely,7)=num((i-1)*nely+1:i*nely,5)+nely;
        num((i-1)*nely+1:i*nely,8)=ydot*(i-1)+xdot*ydot+xdot*nely+1:ydot*(i-1)+xdot*ydot+xdot*nely+nely;
        num((i-1)*nely+1:i*nely,6)=num((i-1)*nely+1:i*nely,8)+1;        
        
    end
    
%% preparation for FEA

nodeloc=zeros(nodenum,2);
for i=1:xdot
    nodeloc(ydot*(i-1)+1:ydot*i,1)=i-1;
    a=ydot-1:-1:0;
    nodeloc(ydot*(i-1)+1:ydot*i,2)=a;
    
    nodeloc(xdot*ydot+nely*(i-1)+1:xdot*ydot+nely*i,1)=i-1;
    b=ydot-1-0.5:-1:0.5;
    nodeloc(xdot*ydot+nely*(i-1)+1:xdot*ydot+nely*i,2)=b;
    
end


for i=1:nelx
    nodeloc(xdot*ydot+nely*xdot+ydot*(i-1)+1:xdot*ydot+nely*xdot+ydot*i,1)=i-1+0.5;
    a=ydot-1:-1:0;
    nodeloc(xdot*ydot+nely*xdot+ydot*(i-1)+1:xdot*ydot+nely*xdot+ydot*i,2)=a;
end



K=zeros(2*nodenum,2*nodenum);
u=zeros(2*nodenum,1);
f=zeros(2*nodenum,1);


%% boundary condition

f(2*xdot*ydot-ydot,1)=F; %%%施加负载
fixednode=union([1:ydot*2],[2*(xdot*ydot+1)-1:2*(xdot*ydot+nely)]); %%%固定约束
allnode=[1:nodenum*2];  
freenode=setdiff(allnode,fixednode);

%% stiffness matrix

for i=1:elemnum
    elemloc=zeros(4,2);  
    elemloc(1:4,:)=nodeloc(num(i,1:4),:);
    ke=elemstiffness(E,nu,t,elemloc);
    K=stiffness(K,ke,num(i,:));
end

u(freenode)=K(freenode,freenode)\f(freenode);

u(fixednode)=0;


%% display result
nodexlocold=nodeloc(1:ydot*xdot,1);
nodeylocold=nodeloc(1:ydot*xdot,2);

nodexloc=nodeloc(1:ydot*xdot,1)+u(1:2:2*ydot*xdot-1);
nodeyloc=nodeloc(1:ydot*xdot,2)+u(2:2:2*ydot*xdot);
elem=num(:,1:4);
display_2D(elem,nodexloc,nodeyloc,nodexlocold,nodeylocold)  

end