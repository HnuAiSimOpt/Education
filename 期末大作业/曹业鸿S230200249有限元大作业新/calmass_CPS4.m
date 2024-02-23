function [Mass] = calmass_CPS4(gcoords,nodes,Density,Thickness,nnpe,nds)

    nes = length(nodes(:,1));
    nns = length(gcoords(:,1));

mass=zeros(nns,nns);
for count=1:nes
   x1=gcoords(nodes(count,2),2);y1=gcoords(nodes(count,2),3);
   x2=gcoords(nodes(count,3),2);y2=gcoords(nodes(count,3),3);
   x3=gcoords(nodes(count,4),2);y3=gcoords(nodes(count,4),3);
   x4=gcoords(nodes(count,5),2);y4=gcoords(nodes(count,5),3);
   area1=abs(0.5*det([1 x1 y1;1 x2 y2;1 x3 y3]));
   area2=abs(0.5*det([1 x1 y1;1 x4 y4;1 x3 y3]));
   area=area1+area2;
   for count2=1:nnpe
       mass(nodes(count,count2),nodes(count,count2))=mass(nodes(count,count2),nodes(count,count2))+area*Density*Thickness/4;
   end
end 
Mass=zeros(nds,nds);
for count=1:nns
    Mass(2*count-1,2*count-1)=mass(count,count);
    Mass(2*count,2*count)=mass(count,count);
end
return