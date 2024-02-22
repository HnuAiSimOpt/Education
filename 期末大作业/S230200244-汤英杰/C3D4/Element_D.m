function D = Element_D(iopt,elastic,poisson)
%----------------------------------------------%
% Purpose:
% determine the constitutive equation for isotropic material
%
% Synopis;
% [matmtrx]=fematiso(iopt,elastic,poison)
%
% Variable Description
% elastic- elastic modulus
% iopt=1-plane stress analysis
% iopt=2-plane strain analysis
% iopt=3-axisymmetric analysis
% iopt=4-three dimensional analysis
%----------------------------------------------%
if iopt==1   %plane stress
    D=elastic/(1-poisson*poisson)*...
        [1 poisson 0;...
        poisson 1 0;...
        0 0 (1-poisson)/2];
    %
elseif iopt==2 %plane strain
    D=elastic/((1+poisson)*(1-2*poisson))*...
        [(1-poisson) poisson 0;
        poisson (1-poisson) 0;
        0 0 (1-2*poisson)/2];
    %
elseif iopt==3 %axisymmetry
    D=elastic/((1+poisson)*(1-2*poisson))*...
        [(1-poisson) poisson poisson 0;
        poisson (1-poisson) poisson 0;
        poisson poisson (1-poisson) 0;
        0 0 0 (1-2*poisson)/2];
    %
else
    D=elastic/((1+poisson)*(1-2*poisson))*...
        [(1-poisson) poisson poisson 0 0 0;
        poisson (1-poisson) poisson 0 0 0;
        poisson poisson (1-poisson) 0 0 0;
        0 0 0 (1-2*poisson)/2 0 0;
        0 0 0 0 (1-2*poisson)/2 0;
        0 0 0 0 0 (1-2*poisson)/2];
    %

end