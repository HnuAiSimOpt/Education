function [K_c, F_c] = Boundary_Apply(K, F, Dof,Constrain_nodes)
%----------------------------------------------%
%function:
%   Apply Boundary conditions on Stiffness K and Load Matrix F
%input:
%   K:
%       System Stiffness Matrix
%   F:
%       System Load Matrix
%   Dof:
%       The Dof should set to be casted 
%   Constrain_nodes:
%       The Dofs set to be zero
%----------------------------------------------%
Constrain_dofs = zeros(length(Constrain_nodes),3);
for i = 1:Dof
    Constrain_dofs(:,i) = (Constrain_nodes-1)*Dof+i;
end
switch Dof
    case 1
        Constrain = Constrain_dofs(:,1);
    case 2
        Constrain = [Constrain_dofs(:,1);Constrain_dofs(:,2)];
    case 3
        Constrain = [Constrain_dofs(:,1);Constrain_dofs(:,2);Constrain_dofs(:,3)];
end
K(Constrain,:) = [];
F(Constrain,:) = [];
K(:,Constrain) = [];
K_c = K;
F_c = F;

end