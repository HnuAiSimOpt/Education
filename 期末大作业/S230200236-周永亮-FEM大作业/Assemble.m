<<<<<<< HEAD
function[kk]=Assemble(kk,nd,k,nnel,ndof)
%组装刚度矩阵
for i=1:nnel
    for j=1:nnel
        kk(ndof*(nd(i)-1)+1:ndof*nd(i),ndof*(nd(j)-1)+1:ndof*nd(j))=...
            kk(ndof*(nd(i)-1)+1:ndof*nd(i),ndof*(nd(j)-1)+1:ndof*nd(j))+...
            k(ndof*(i-1)+1:ndof*i,ndof*(j-1)+1:ndof*j);%将单元刚度矩阵放到总体刚度矩阵中
    end
=======
function[kk]=Assemble(kk,nd,k,nnel,ndof)
%组装刚度矩阵
for i=1:nnel
    for j=1:nnel
        kk(ndof*(nd(i)-1)+1:ndof*nd(i),ndof*(nd(j)-1)+1:ndof*nd(j))=...
            kk(ndof*(nd(i)-1)+1:ndof*nd(i),ndof*(nd(j)-1)+1:ndof*nd(j))+...
            k(ndof*(i-1)+1:ndof*i,ndof*(j-1)+1:ndof*j);%将单元刚度矩阵放到总体刚度矩阵中
    end
>>>>>>> 0c2441c35348cc2cb0d4bfa4f93dc736c7d002fd
end