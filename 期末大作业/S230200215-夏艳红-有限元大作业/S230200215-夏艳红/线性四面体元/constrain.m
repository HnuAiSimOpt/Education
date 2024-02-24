function cons=constrain(Boundary_nodes,ndof)
con_dofs = zeros(length(Boundary_nodes),3);
for i = 1:ndof
    con_dofs(:,i) = (Boundary_nodes-1)*ndof+i;
end
switch ndof
    case 1
        cons = con_dofs(:,1);
    case 2
        cons = [con_dofs(:,1);con_dofs(:,2)];
    case 3
        cons = [con_dofs(:,1);con_dofs(:,2);con_dofs(:,3)];
end