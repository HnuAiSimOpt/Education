<<<<<<< HEAD
function[kk]=Boundary(bcdof,nbc,kk,ndof)
%添加边界条件
for i=1:nbc
    for j=ndof*(bcdof(i)-1)+1:ndof*bcdof(i)
        kk(j,:)=0;
        kk(:,j)=0;
        kk(j,j)=1;
    end
=======
function[kk]=Boundary(bcdof,nbc,kk,ndof)
%添加边界条件
for i=1:nbc
    for j=ndof*(bcdof(i)-1)+1:ndof*bcdof(i)
        kk(j,:)=0;
        kk(:,j)=0;
        kk(j,j)=1;
    end
>>>>>>> 0c2441c35348cc2cb0d4bfa4f93dc736c7d002fd
end