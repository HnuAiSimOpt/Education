function [matmtrx]=fematiso(iopt,elastic,poisson)
if  iopt==1        % plane stress
 matmtrx= elastic/(1-poisson*poisson)* ...
    [1  poisson 0; ...
     poisson  1  0; ...
     0  0  (1-poisson)/2];
elseif   iopt==2        % plane strain
 matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson 0; 
   poisson  (1-poisson)  0;
   0  0  (1-2*poisson)/2];
 elseif  iopt==3     % axisymmetry
 matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson  0; 
   poisson  (1-poisson)   poisson  0;
   poisson  poisson  (1-poisson)   0;
   0    0    0   (1-2*poisson)/2];
elseif  iopt==4        % three-dimension
 matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson   0   0    0; 
   poisson  (1-poisson)   poisson   0   0    0;
   poisson  poisson  (1-poisson)    0   0    0;
   0    0    0    (1-2*poisson)/2   0    0;
   0    0    0    0    (1-2*poisson)/2   0;
   0    0    0    0    0   (1-2*poisson)/2];
end
return
