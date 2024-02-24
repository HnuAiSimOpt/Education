function [shape,nderiv] = shapeFunct_Beam(xi)
%%%% shapeN : 形函数 N1 和 N2
%%%% N_deriv: 推导 N1 and N2 w.r.t. xi 
%%%% xi: 自然坐标 (-1 ... +1)

shape =1/2*[1-xi;1+xi];

nderiv = [-1;1]/2;

end
  
