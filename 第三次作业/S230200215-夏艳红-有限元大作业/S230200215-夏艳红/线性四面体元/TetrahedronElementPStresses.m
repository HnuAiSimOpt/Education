function y=TetrahedronElementPStresses(sigma)
s1=sigma(1)+sigma(2)+sigma(5);
s2=sigma(1)*sigma(2)+sigma(1)*sigma(3)+sigma(2)*sigma(3)-sigma(4)*sigma(4)-sigma(5)*sigma(5)-sigma(6)*sigma(6);
ms3=[sigma(1) sigma(4) sigma(6);sigma(4) sigma(2) sigma(5);sigma(6) sigma(5) sigma(3)];
s3=det(ms3);
y=[s1;s2;s3];