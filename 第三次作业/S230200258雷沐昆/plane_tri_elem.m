function Ke=plane_tri_elem(D,h,coords)

b=zeros(1,3);c=zeros(1,3);
coord=[coords;coords];

for k=1:3
    b(k)=coord(k+1,2)-coord(k+2,2);
    c(k)=-coord(k+1,1)+coord(k+2,1);
end

area=0.5*det([ones(3,1),coords]);

B=zeros(3,6);
B(1,1:2:end)=b;
B(2,2:2:end)=c;
B(3,1:2:end)=c;
B(3,2:2:end)=b;
B=B/(2*area);

Ke=B'*D*B*h*area;

