function [Kglobal,Fglobal]=plane_assembly_tri(Kglobal,Fglobal,node,element,bodyforce,D,h)

numelem=size(element,1);
for iel=1:numelem
    nod=element(iel,:);
    sctr=zeros(1,6);
    sctr(1:2:end)=2*nod-1;
    sctr(2:2:end)=2*nod  ;
    coords=node(nod,:);
    center=mean(coords,1);
    area=0.5*det([ones(3,1),coords]);
    
    Ke=plane_tri_elem(D(center),h,coords);
    Fe=[bodyforce(center);bodyforce(center);bodyforce(center)]*area*h/3;
    
    Kglobal(sctr,sctr)=Kglobal(sctr,sctr)+Ke;
    Fglobal(sctr)     =Fglobal(sctr)+Fe;
end
    
    
