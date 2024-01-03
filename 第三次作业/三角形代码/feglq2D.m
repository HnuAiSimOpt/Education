function [pointArray,weightArray]=feglq2D(nglx,ngly)
if nglx>ngly
    ngl=nglx;
else
    ngl=ngly;
end
pointArray = zeros(ngl,2);
weightArray = zeros(ngl,2);
[pointx,weightx] = feglq1D(nglx);
[pointy,weighty] = feglq1D(ngly); 
for countx=1:nglx
     pointArray(countx,1) = pointx(countx);
    weightArray(countx,1) = weightx(countx);
end
for county=1:ngly
     pointArray(county,2) = pointy(county);
    weightArray(county,2) = weighty(county);
end
return