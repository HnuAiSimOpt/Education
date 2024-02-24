function[kk]=boundary(bcdof,nbc,kk,ndof)
% 施加边界条件
for i=1:nbc
    for j=ndof*(bcdof(i)-1)+1:ndof*bcdof(i)
        kk(j,:)=0;
        kk(:,j)=0;
        kk(j,j)=1;
    end
end