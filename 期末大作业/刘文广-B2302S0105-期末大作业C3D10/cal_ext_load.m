function [ff]=cal_ext_load(gcoord,nodes,ndof,L,p0)
%--------------------------------------------------------------------------
% Purpose:
%     Calculate system force vector
%
% Synopsis:
%     [ff]=cal_ext_load(gcoord,nodes,ndof,L,p0)
%
% Variable descriptions:
%     ff - system force vector
%     gcoord - nodal coordinate values for each node 
%     nodes - nodal connectivity for each element
%     ndof - number of dofs per node 
%     L - length of the cantilever beam in x-axis 
%     p0 - magnitude of the external load 
%--------------------------------------------------------------------------
nnode=length(gcoord);
[nel,nnel]=size(nodes);
sdof=nnode*ndof;
edof=nnel*ndof;
ff=sparse(sdof,1);
erro=1.0e-4;

for iel=1:nel
    k=0;
    f=zeros(edof,1);
    nd=zeros(nnel,1);
    xcoord=zeros(nnel,1);
    ycoord=zeros(nnel,1);
    zcoord=zeros(nnel,1);
    mat=zeros(6,4);
    for i=1:nnel
        nd(i)=nodes(iel,i);
        xcoord(i)=gcoord(nd(i),1);
        ycoord(i)=gcoord(nd(i),2);
        zcoord(i)=gcoord(nd(i),3);
        if (abs(xcoord(i)-L)<=erro)
            k=k+1;
            mat(k,1)=i;
            mat(k,2)=xcoord(i);
            mat(k,3)=ycoord(i);
            mat(k,4)=zcoord(i);
        end
    end
    
    if k==6
        x=[mat(1,2);mat(2,2);mat(3,2)];
        y=[mat(1,3);mat(2,3);mat(3,3)];
        z=[mat(1,4);mat(2,4);mat(3,4)];
        
        [area]=cal_area_3T(x,y,z);
        fe=zeros(3*ndof,1);
        fe(2)=area*p0/3;
        fe(5)=area*p0/3;
        fe(8)=area*p0/3;
        [index]=feeldof(mat(1:3,1),3,ndof);
        [f]=feasmbl_f(f,fe,index);
    end
    [index]=feeldof(nd,nnel,ndof);
    [ff]=feasmbl_f(ff,f,index);
end

return