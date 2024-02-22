function [kk]=cal_stiff_t10(gcoord,nodes,ndof,matmtx)
%--------------------------------------------------------------------------
% Purpose:
%     Calculate system stiffness matrix using FEM-T10 element
%
% Synopsis:
%     kk - system stiffness matrix 
%     gcoord - nodal coordinate values for each node 
%     nodes - nodal connectivity for each element
%     ndof - number of dofs per node 
%     matmtx - the constitutive matrix for 3D-solid
%--------------------------------------------------------------------------
nnode=length(gcoord);
[nel,nnel]=size(nodes);
edof=nnel*ndof;
sdof=nnode*ndof;
kk=sparse(sdof,sdof);

iham=4;                                  % number of Hammer integration points(=1;=4;=5)
[points,weights]=Hammer_T4(iham);        % points and weights

for iel=1:nel
    ke=zeros(edof,edof);
    nd=zeros(nnel,1);
    xcoord=zeros(nnel,1);
    ycoord=zeros(nnel,1);
    zcoord=zeros(nnel,1);
    for i=1:nnel
        nd(i)=nodes(iel,i);
        xcoord(i)=gcoord(nd(i),1);
        ycoord(i)=gcoord(nd(i),2);
        zcoord(i)=gcoord(nd(i),3);
    end
    x=[xcoord(1);xcoord(2);xcoord(3);xcoord(4)];
    y=[ycoord(1);ycoord(2);ycoord(3);ycoord(4)];
    z=[zcoord(1);zcoord(2);zcoord(3);zcoord(4)];
    [a,b,c,d,vol]=cal_der_t10(x,y,z);
    
    for ih=1:iham
        L=zeros(1,nnel);
        L(1)=points(ih,1);
        L(2)=points(ih,2);
        L(3)=points(ih,3);
        L(4)=points(ih,4);
        wt=weights(ih);
        
        [shap,bmat]=cal_shap_bmat_t10(nnel,ndof,L,a,b,c,d,vol);
        ke=ke+6*vol*bmat'*matmtx*bmat*wt/6;
    end
    [index]=feeldof(nd,nnel,ndof);
    [kk]=feasmbl_k(kk,ke,index);
end

return
        
        
        





