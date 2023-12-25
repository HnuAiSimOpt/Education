function [index] = feeldof(node,nnpe,ndpn)
index = [];
for count = 1:nnpe
    for countin = 1:ndpn
        index = [index ndpn*node(count)-(ndpn-countin)];
    end
end
end