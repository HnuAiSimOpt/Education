
function [Ke]=Ke_C3D4(D, ElementNodeCoordinate)
%初始化单元刚度阵
Ke=zeros(12,12);
            % 计算形函数对总体坐标的导数（NDerivative）及雅可比矩阵行列式（JacobiDET）
            [NDxyz, JacobiDET] = ShapeFunction( ElementNodeCoordinate);%[DN1Dx DN2Dx DN3Dx;DN1Dy DN2Dy DN3Dy;……]
            Ve = JacobiDET/6;%
            %计算B矩阵  利用形函数对总体坐标的导数（NDxyz）对B进行计算
            B=zeros(6,12);
            for i=1:4
                sub=(i-1)*3+1:(i-1)*3+3;%子矩阵范围
                B(:,sub)=[NDxyz(1,i) 0         0;%NDx
                    0         NDxyz(2,i) 0;%NDy
                    0         0         NDxyz(3,i);%NDz
                    NDxyz(2,i) NDxyz(1,i) 0;
                    0         NDxyz(3,i) NDxyz(2,i);
                    NDxyz(3,i) 0         NDxyz(1,i)];
            end
            Ke=Ke+Ve*B'*D*B;  %数值积分，多个积分点叠加
end