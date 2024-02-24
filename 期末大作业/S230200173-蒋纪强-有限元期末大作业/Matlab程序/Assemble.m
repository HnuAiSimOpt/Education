function[kk]=Assemble(kk,nd,k,nnel,ndof)
%×é×°¸Õ¶È¾ØÕó
for i=1:nnel
    for j=1:nnel
        kk(ndof*(nd(i)-1)+1:ndof*nd(i),ndof*(nd(j)-1)+1:ndof*nd(j))=...
            kk(ndof*(nd(i)-1)+1:ndof*nd(i),ndof*(nd(j)-1)+1:ndof*nd(j))+...
            k(ndof*(i-1)+1:ndof*i,ndof*(j-1)+1:ndof*j);
    end
end