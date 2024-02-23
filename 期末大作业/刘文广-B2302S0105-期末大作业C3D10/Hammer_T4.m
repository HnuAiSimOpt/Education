function [points,weights]=Hammer_T4(iham)
%--------------------------------------------------------------------------
% Purpose:
%     Input the integration points and weights for Hammer integration
%
% Synopsis:
%     [points,weights]=Hammer_T4(iham)
%
% Variable descriptions:
%     points - coordinate values of integration points 
%     weights - weights of integration points 
%     iham - number of integration points 
%--------------------------------------------------------------------------

if iham==1
    points=[0.25 0.25 0.25 0.25];
    weights=1;
elseif iham==4
    a=0.58541020;
    b=0.13819660;
    points=[a b b b;b a b b;b b a b;b b b a];
    weights=[0.25;0.25;0.25;0.25];
elseif iham==5
    points=[1/4 1/4 1/4 1/4;...
            1/2 1/6 1/6 1/6;...
            1/6 1/2 1/6 1/6;...
            1/6 1/6 1/2 1/6;...
            1/6 1/6 1/6 1/2];
    weights=[-4/5;9/20;9/20;9/20;9/20];
end

return