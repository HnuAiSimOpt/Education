
function K=stiffness(K,ke,elemnode)

Knum=zeros(1,6);
Knum(1:2:5)=2*elemnode(:)-1;
Knum(2:2:6)=2*elemnode(:);

for i=1:6
   for j=1:6
       K(Knum(i),Knum(j))=K(Knum(i),Knum(j))+ke(i,j);
   end
end
end
