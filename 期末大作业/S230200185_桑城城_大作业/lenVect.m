function [lenVect] = lenVect(Vect)
%%%% nPointGen = 需要生成的材料点数
lenVect = (Vect(1,1)^2+Vect(1,2)^2+Vect(1,3)^2)^(1/2);
end

