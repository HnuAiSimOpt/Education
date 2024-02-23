function [detJ, invJ] = detJ_invJ(xyz,eNodei)
%%% 计算单元长度, det(J) and J^-1
Le = sqrt((xyz(eNodei(2),1)-xyz(eNodei(1),1))^2 + ...
        (xyz(eNodei(2),2)-xyz(eNodei(1),2))^2 + ...
        (xyz(eNodei(2),3)-xyz(eNodei(1),3))^2);

    %%%% det(J) and J^-1    
detJ = Le/2;
invJ = 1/detJ;

end

