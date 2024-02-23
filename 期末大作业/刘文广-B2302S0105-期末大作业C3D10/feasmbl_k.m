function [kk]=feasmbl_k(kk,k,index)
%--------------------------------------------------------------------------
% Purpose:
%     Assembly of element stiffness matrix into system stiffness matrix 
%
% Synopsis:
%     [kk]=feasmbl_k(kk,k,index)
%
% Variable descriptions:
%     kk - system stiffness matrix 
%     k - element stiffness matrix 
%     index - d.o.f. vector associated with an element 
%--------------------------------------------------------------------------

edof=length(index);

for i=1:edof
    ii=index(i);
    for j=1:edof
        jj=index(j);
        kk(ii,jj)=kk(ii,jj)+k(i,j);
    end
end

return