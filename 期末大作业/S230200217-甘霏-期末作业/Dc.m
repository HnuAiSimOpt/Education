function [D]=Dc(iopt,elastic,poisson)
% µ√µΩµØ–‘æÿ’ÛD
if iopt==1
    D=elastic/(1-poisson^2)*...
        [1 poisson 0
        poisson 1 0
        0 0 (1-poisson)/2];
elseif iopt==2
    D=elastic/((1+poisson)*(1-2*poisson))*...
        [1-poisson poisson 0
        poisson 1-poisson 0
        0 0 (1-2*poisson)/2];
elseif iopt==3
    D=elastic/((1+poisson)*(1-2*poisson))*...
        [1-poisson poisson poisson 0
        poisson 1-poisson poisson 0
        poisson poisson 1-poisson 0
        0 0 0 (1-2*poisson)/2];
else
    D=elastic/((1+poisson)*(1-2*poisson))*...
        [1-poisson poisson poisson 0 0 0
        poisson 1-poisson poisson 0 0 0
        poisson poisson 1-poisson 0 0 0
        0 0 0 (1-2*poisson)/2 0 0
        0 0 0 0 (1-2*poisson)/2 0
        0 0 0 0 0 (1-2*poisson)/2];
end