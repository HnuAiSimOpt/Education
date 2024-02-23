function K=stiffness(K,ke,elemnode)

Knum=zeros(1,16);
Knum(1:2:15)=2*elemnode(:)-1;
Knum(2:2:16)=2*elemnode(:);

for i=1:16
   for j=1:16
       K(Knum(i),Knum(j))=K(Knum(i),Knum(j))+ke(i,j);
   end
end
end