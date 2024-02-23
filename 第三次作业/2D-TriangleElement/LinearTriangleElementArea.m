function y = LinearTriangleElementArea(xi, yi, xj, yj, xm, ym)
%LinearTriangleElementArea      函数用于求解(xi, yi), (xj, yj), (xm, ym)
%                               节点组成的三角形单元面积

y = (xi*(yj - ym) + xj*(ym - yi) + xm*(yi - yj))/2;

end