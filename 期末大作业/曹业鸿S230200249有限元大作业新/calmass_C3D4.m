function [Mass] = calmass_C3D4(gcoords,nodes,Density,Thickness,nnpe,nds)
ndpn=3;
nes=length(nodes(:,1));
Mass=sparse(nds,nds);
x=zeros(nnpe,1);
y=zeros(nnpe,1);
z=zeros(nnpe,1);
for count=1:nes
    cnode=nodes(count,2:nnpe+1);
    for count2=1:nnpe
        x(count2)=gcoords(cnode(count2),2);
        y(count2)=gcoords(cnode(count2),3);
        z(count2)=gcoords(cnode(count2),4);
    end
    V6=[x(2)-x(1),y(2)-y(1),z(2)-z(1);
        x(3)-x(1),y(3)-y(1),z(3)-z(1);
        x(4)-x(1),y(4)-y(1),z(4)-z(1);];
    V=det(V6)/6;
    masse=abs(V)*Density;
    Masse=0.25*eye(nnpe*ndpn)*masse;
    index=feeldof(cnode,nnpe,ndpn);
    Mass = assemblekk(Mass,Masse,index);
end
return