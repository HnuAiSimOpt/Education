function [bcdof,bcval]=In_bcdof(bound,ndof)
%--------------------------------------------------------------------------
% Purpose:
%     Input data for boundary conditions
%
% Synopsis:
%     [bcdof,bcval]=In_bcdof(bound,ndof)
%
% Variable descriptions:
%     bcdof - a vector containing the constrained dofs 
%     bcval - a vector containing the constrained values 
%     bound - boundary node and its coordinates 
%     ndof - number of dofs per node 
%--------------------------------------------------------------------------
[nbcnode,dim]=size(bound);
bcdof=zeros(1,nbcnode*ndof);

for inode=1:nbcnode
    node=bound(inode,1);
    bcdof(inode*ndof-2)=node*ndof-2;
    bcdof(inode*ndof-1)=node*ndof-1;
    bcdof(inode*ndof-0)=node*ndof-0;
end

nbcdof=length(bcdof);
bcval=zeros(1,nbcdof);

return
