function [kk]= assemblekk(kk,k,index)

Len = length(index); 
for n1=1:Len
    for n2=1:Len
        kk(index(n1),index(n2))= kk(index(n1),index(n2))+k(n1,n2);
    end
end