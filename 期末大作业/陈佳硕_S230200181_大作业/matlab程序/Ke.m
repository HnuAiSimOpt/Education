%%%%һ�������嵥Ԫ��λ�նȾ���
% Ke��λ�նȾ���
% D�������ͬ���ߵ��Բ���Ӧ��-Ӧ�����  
% ElementNodeCoordinate��Ԫ�ڵ����꣨8*3��ÿһ�д���һ���ڵ�����꣩
function [Ke]=Ke(D, ElementNodeCoordinate)
% ��˹���ֵ�����
GaussCoordinate=[-0.57735026918963D0, 0.57735026918963D0];
%��˹���ֵ�Ȩ��
GaussWeight=[1.00000000000000D0, 1.00000000000000D0];
%��ʼ����Ԫ�ն���
Ke=zeros(24,24);
%ѭ����˹��
for X=1:2
    for Y=1:2
        for Z=1:2
            GP1=GaussCoordinate(X); GP2=GaussCoordinate(Y); GP3=GaussCoordinate(Z); %��˹������
            % �����κ�������������ĵ�����NDerivative�����ſɱȾ�������ʽ��JacobiDET��
            [~,NDerivative, JacobiDET] = ShapeFunction([GP1 GP2 GP3], ElementNodeCoordinate);
            Coefficient=GaussWeight(X)*GaussWeight(Y)*GaussWeight(Z)*JacobiDET;
            %����B����  �����κ�������������ĵ�����NDerivative����B���м���
            B=zeros(6,24);
            for I=1:8
                COL=(I-1)*3+1:(I-1)*3+3;
                B(:,COL)=[NDerivative(1,I) 0         0;
                    0         NDerivative(2,I) 0;
                    0         0         NDerivative(3,I);
                    NDerivative(2,I) NDerivative(1,I) 0;
                    0         NDerivative(3,I) NDerivative(2,I);
                    NDerivative(3,I) 0         NDerivative(1,I)];
            end
            Ke=Ke+Coefficient*B'*D*B;  %���Ӹն���
        end
    end
end
end