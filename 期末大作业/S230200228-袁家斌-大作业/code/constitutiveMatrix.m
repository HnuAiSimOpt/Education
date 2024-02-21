function [constitutive] = constitutiveMatrix(E, mu)

%% 返回材料的面内本构矩阵
constitutive.in = E/(1-mu^2)*[1, mu, 0
                            mu, 1, 0
                            0 0 0.5*(1-mu)];

%% 返回材料的面外剪切本构矩阵
G = E /2 /(1 + mu);
constitutive.out = G * eye(2);

%% 厚板修正系数
constitutive.alpha = 5.0 / 6.0;
end

