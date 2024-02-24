function [kk,ff] = loadbc(kk,BC,ff)
for count = 1:length(BC)
    x = BC(count);
    for count_in = 1:length(kk(1,:))
        kk(x,count_in) = 0;
        kk(count_in,x) = 0;
    end
    kk(x,x)=1;
    ff(x)=0;
end
end