
function [Ke]=Ke_C3D4(D, ElementNodeCoordinate)
%��ʼ����Ԫ�ն���
Ke=zeros(12,12);
            % �����κ�������������ĵ�����NDerivative�����ſɱȾ�������ʽ��JacobiDET��
            [NDxyz, JacobiDET] = ShapeFunction( ElementNodeCoordinate);%[DN1Dx DN2Dx DN3Dx;DN1Dy DN2Dy DN3Dy;����]
            Ve = JacobiDET/6;%
            %����B����  �����κ�������������ĵ�����NDxyz����B���м���
            B=zeros(6,12);
            for i=1:4
                sub=(i-1)*3+1:(i-1)*3+3;%�Ӿ���Χ
                B(:,sub)=[NDxyz(1,i) 0         0;%NDx
                    0         NDxyz(2,i) 0;%NDy
                    0         0         NDxyz(3,i);%NDz
                    NDxyz(2,i) NDxyz(1,i) 0;
                    0         NDxyz(3,i) NDxyz(2,i);
                    NDxyz(3,i) 0         NDxyz(1,i)];
            end
            Ke=Ke+Ve*B'*D*B;  %��ֵ���֣�������ֵ����
end