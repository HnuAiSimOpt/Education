

function [NDxyz,JacobiDET] = ShapeFunction(ElementNodeCoordinate)
%计算形函数及形函数对局部坐标ksi eta zeta的导数
NDL = [-1 1 0 0;-1 0 1 0;-1 0 0 1];%3*4  [N1Dksi N2Dksi N3Dksi N4Dksi;N1Deta N2Deta N3Deta N4Deta……]
Jacobi = NDL*ElementNodeCoordinate;%计算雅可比矩阵3*4 4*3
JacobiDET = det(Jacobi);%计算雅可比行列式3*3 [DxDksi DyDksi DzDksi;DxDeta……
JacobiINV=inv(Jacobi);%对雅可比行列式求逆3*3
NDxyz=JacobiINV*NDL;%利用雅可比行列式的逆计算形函数对结构坐标的导数[DN1Dx DN2Dx DN3Dx;DN1Dy DN2Dy DN3Dy;……]
end