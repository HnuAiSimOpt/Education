function [point3, weight3] = GLI_PW3(nglx, ngly, nglz)

if nglx>ngly
    if nglx>ngz
        ngl = nglx;
    else
        ngl = nglz;
    end
else
    if ngly >nglz
        ngl = ngly;
    else
        ngl = nglz;
    end
end

point3 = zeros(ngl,3);
weight3 = zeros(ngl,3);

[pointx, weightx] = GLI_PW1(nglx);
[pointy, weighty] = GLI_PW1(ngly);
[pointz, weightz] = GLI_PW1(nglz);

for intx = 1:nglx
    point3(intx, 1) = pointx(intx);
    weight3(intx, 1) = weightx(intx);
end
for inty = 1:ngly
    point3(inty, 2) = pointy(inty);
    weight3(inty, 2) = weighty(inty);
end
for intz = 1:nglz
    point3(intz, 3) = pointz(intz);
    weight3(intz, 3) = weightz(intz);
end

end