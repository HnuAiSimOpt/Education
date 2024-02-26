function D = MatrixD(E,nu)
%计算弹性矩阵
D = E/(1+nu)/(1-2*nu)*[1-nu, nu, nu, 0, 0, 0;
                                       nu, 1-nu, nu, 0, 0, 0;
                                       nu, nu, 1-nu, 0, 0, 0;
                                       0, 0, 0, (1-2*nu)/2, 0, 0;
                                       0, 0, 0, 0, (1-2*nu)/2, 0;
                                       0, 0, 0, 0, 0, (1-2*nu)/2];
end