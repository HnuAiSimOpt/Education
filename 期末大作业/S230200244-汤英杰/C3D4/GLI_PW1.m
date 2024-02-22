function [point1, weight1] = GLI_PW1(ngl)
%----------------------------------------------%
%function:
%   return Gauss-Legendre quadrature for intergation
%input:
%   ngl:
%       number of intergation points
%       from 1 to 5
%output:
%   point1:
%       points of intergation point
%   weights:
%       weights of intergation points coordiante to point1
%----------------------------------------------%
point1 = zeros(ngl, 1);
weight1 = zeros(ngl,1);
switch(ngl)
    case 1
        point1 = 0;
        weight1 = 2;
    case 2
        point1 = [-0.577350269189626;0.577350269189626];
        weight1 = [1;1];
    case 3
        point1 = [-0.774596669241483;0;0.774596669241483];
        weight1 = [0.55555555;0.88888888;0.55555555];
    case 4
        point1 = [-0.861136311594053;-0.339981043584856;0.339981043584856;0.861136311594053];
        weight1 = [0.347854845137454;0.652145154862546;0.652145154862546;0.347854845137454];   
    case 5
        point1 = [-0.906179845938664;-0.538469310105683;0;0.538469310105683;0.906179845938664];
        weight1 = [0.236926885056189;0.478628670499366;0.568888888888889;0.478628670499366;0.236926885056189];   
end

end