function [PointArray,WeightArray] = feglq1D(ngl)
 PointArray = zeros(ngl,1);
WeightArray = zeros(ngl,1);
if ngl==1
     PointArray(1) = 0.0;
    WeightArray(1) = 2.0;
elseif ngl==2 
     PointArray(1) = -0.577350269189626;
     PointArray(2) = -PointArray(1);
    WeightArray(1) = 1.0;
    WeightArray(2) = WeightArray(1);
elseif ngl==3
     PointArray(1) = -0.774596669241483;
     PointArray(2) = 0.0;
     PointArray(3) = -PointArray(1);
    WeightArray(1) = 0.555555555555556;
    WeightArray(2) = 0.888888888888889;
    WeightArray(3) = WeightArray(1);
elseif ngl==4
     PointArray(1) = -0.861136311594053;
     PointArray(2) = -0.339981043584856;
     PointArray(3) = -PointArray(2);
     PointArray(4) = -PointArray(1);
    WeightArray(1) = 0.347854845137454;
    WeightArray(2) = 0.652145154862546;
    WeightArray(3) = WeightArray(2);
    WeightArray(4) = WeightArray(1);
elseif ngl==5
     PointArray(1) = -0.906179845938664;
     PointArray(2) = -0.538469310105683;
     PointArray(3) = 0.0;
     PointArray(4) = -PointArray(2);
     PointArray(5) = -PointArray(1);
    WeightArray(1) = 0.236926885056189;
    WeightArray(2) = 0.478628670499366;
    WeightArray(3) = 0.568888888888889;
    WeightArray(4) = WeightArray(2);
    WeightArray(5) = WeightArray(1);
else
    fprintf("Err: Not Supported!\n");
    fprintf("Number of Gauss-Legendre - NO MORE THAN 5!\n");
end
return