function disp = solution(nDof,fixedDof,K,force)
activeDof = setdiff((1:nDof)',fixedDof);
U = K(activeDof,activeDof)\force(activeDof);
disp = zeros(nDof,1);
disp(activeDof) = U;