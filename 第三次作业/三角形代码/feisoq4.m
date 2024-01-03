function [shapeq4,dhdxiq4,dhdetaq4]=feisoq4(xi,eta)
shapeq4(1) = 0.25*(1-xi)*(1-eta);
shapeq4(2) = 0.25*(1+xi)*(1-eta);
shapeq4(3) = 0.25*(1+xi)*(1+eta);
shapeq4(4) = 0.25*(1-xi)*(1+eta);
dhdxiq4(1) = -0.25*(1-xi);
dhdxiq4(2) = 0.25*(1-xi);
dhdxiq4(3) = 0.25*(1+xi);
dhdxiq4(4) = -0.25*(1+xi);
dhdetaq4(1) = -0.25*(1-eta);
dhdetaq4(2) = -0.25*(1+eta);
dhdetaq4(3) =  0.25*(1+eta);
dhdetaq4(4) =  0.25*(1-eta);
end
