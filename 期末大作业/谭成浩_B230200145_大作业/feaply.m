function[kk,ff]=feaply(kk,ff,bcdof,bcval)
n=length(bcdof);
sdof=size(kk);
for i=1:n
    c=bcdof(i);
    for j=1:sdof
        kk(c,j)=0;
        kk(j,c)=0;
    end
    kk(c,c)=1;
    ff(c)=bcval(i);
end