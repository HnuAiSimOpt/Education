function [J]=Jacobain(nnel,dndr,dnds,dndt,X,Y,Z)
%º∆À„—≈ø…±»æÿ’Û
J=zeros(3,3);
for i=1:nnel
    J(1,1)=J(1,1)+dndr(i)*X(i);
    J(1,2)=J(1,2)+dndr(i)*Y(i);
    J(1,3)=J(1,3)+dndr(i)*Z(i);
    J(2,1)=J(2,1)+dnds(i)*X(i);
    J(2,2)=J(2,2)+dnds(i)*Y(i);
    J(2,3)=J(2,3)+dnds(i)*Z(i);
    J(3,1)=J(3,1)+dndt(i)*X(i);
    J(3,2)=J(3,2)+dndt(i)*Y(i);
    J(3,3)=J(3,3)+dndt(i)*Z(i);
end