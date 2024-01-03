function y = LinearTriangleElementStiffness(E,NU,t,xi,yi,xj,yj,xm,ym,p)

A=(xi*(yj-ym)+xj*(ym-yi)+xm*(yi-yj))/2;%计算单元面积

betai=yj-ym; %计算单元应变矩阵中的元素

betaj=ym-yi;%计算单元应变矩阵中的元素

betam=yi-yj;%计算单元应变矩阵中的元素

gammai=xm-xj;%计算单元应变矩阵中的元素

gammaj=xi-xm;%计算单元应变矩阵中的元素

gammam=xj-xi;%计算单元应变矩阵中的元素

B=[betai 0 betaj 0 betam 0;%形成单元应变矩阵B

    0 gammai 0 gammaj 0 gammam;

    gammai betai gammaj betaj gammam betam]/(2*A);

if p==1%对于平面应力问题形成弹性矩阵D

    D=(E/(1-NU^2))*[1 NU 0;NU 1 0;0 0 (1-NU)/2];

elseif p==2%对于平面应变问题形成弹性矩阵D

    D=(E/(1+NU)/(1-2*NU))*[1-NU NU 0;NU 1-NU 0;0 0 (1-2*NU)/2];
end

y=t*A*B'*D*B;%根据虚功原理形成单元刚度矩阵