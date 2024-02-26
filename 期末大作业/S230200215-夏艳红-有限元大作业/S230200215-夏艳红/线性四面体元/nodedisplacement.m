function U=nodedisplacement(U_,Boundary_nodes)
U = U_;
for i = 1:length(Boundary_nodes)
    index = Boundary_nodes(i);
    forward_ = U(1:(index-1)*3,:);
    backward_ = U((index-1)*3+1:end,:); 
    U = [forward_;0;0;0;backward_];
end