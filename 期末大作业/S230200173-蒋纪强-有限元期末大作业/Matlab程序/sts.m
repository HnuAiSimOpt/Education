function  z=sts(x)
z=(2/9)*( (x(1)-x(2))^2+(x(2)-x(3))^2+(x(3)-x(1))^2+6*(x(4)^2+x(5)^2+x(6)^2) );
z=sqrt(z); 
end