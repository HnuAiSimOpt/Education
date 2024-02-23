function [bound]=find_node(gcoord,B)
%--------------------------------------------------------------------------
% Purpose:
%     Find the boundary nodes and symmetrical nodes
%
% Synopsis:
%     [bound,sym]=find_node(gcoord,B)
%
% Variable descriptions:
%     bound - boundary node and its coordinate values 
%     sym - symmetrical node and its coordinate values 
%     B - width of the cantilever beam
%     gcoord - nodal coordinate values for each node 
%--------------------------------------------------------------------------
nnode=length(gcoord);
bound=[];
symm=[];
k1=0;
k2=0;
erro=1.0e-3;

for inode=1:nnode 
    x=gcoord(inode,1);
    y=gcoord(inode,2);
    z=gcoord(inode,3);
    if (abs(x-0.0)<=erro)
        k1=k1+1;
        bound(k1,1)=inode;
        bound(k1,2)=x;
        bound(k1,3)=y;
        bound(k1,4)=z;
    end
end
return