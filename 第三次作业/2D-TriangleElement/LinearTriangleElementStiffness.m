function y = LinearTriangleElementStiffness(E, NU, t, xi, yi, xj, yj, xm, ym, P)
%LinearTriangleStiffness    函数用于计算( xi, yi, xj, yj, xm, ym)组成单元的刚度矩阵
%                           E:弹性模量
%                           NU:泊松比
%                           t: 厚度
%                           ( xi, yi, xj, yj, xm, ym):三角形单元坐标
%                           p:选择应变模式（1：平面应力；2：平面应变）
A = LinearTriangleElementArea(xi, yi, xj, yj, xm, ym);
if A<0
   error("A 面积小于0"); 
end
batai = yj - ym;
bataj = ym - yi;
batam = yi - yj;
gamai = xm - xj;
gamaj = xi - xm;
gamam = xj - xi;

B = 1/(2*A)*[
    batai 0 bataj 0 batam 0;
    0 gamai 0 gamaj 0 gamam;
    gamai batai gamaj bataj gamam batam
];
switch P
    case 1
        D = (E/(1-NU^2))*[
            1 NU 0;
            NU 1 0;
            0 0 (1-NU)/2;
        ];
    case 2
        D = (E/((1-2*NU))*(1+NU))*[
            1-NU NU 0;
            NU 1-NU 0;
            0 0 (1-2*NU)/2;
        ];
    otherwise
end

y = t*A*B'*D*B;

end