function F = Load_Apply(Load_nodes, Nodes_num, Dof, Total_Force)
%----------------------------------------------%
%function:
%   Apply Load conditions on Matrix F
%input:
%   Load_nodes:
%       The index of Load nodes
%   Nodes_num:
%       The number of Nodes
%   Dof:
%       The Dof of each Nodes
%   Total_Force:
%       The Total Force applied on Load_nodes
%----------------------------------------------%
F = zeros(Dof*Nodes_num,1);
Load_nodes_num = length(Load_nodes);
for i = 1 :Dof
    F((Load_nodes-1)*Dof+i,1) = Total_Force(i)/Load_nodes_num;

end