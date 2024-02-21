function [constitutive] = constitutiveMatrix(E, mu)

%% ���ز��ϵ����ڱ�������
constitutive.in = E/(1-mu^2)*[1, mu, 0
                            mu, 1, 0
                            0 0 0.5*(1-mu)];

%% ���ز��ϵ�������б�������
G = E /2 /(1 + mu);
constitutive.out = G * eye(2);

%% �������ϵ��
constitutive.alpha = 5.0 / 6.0;
end

