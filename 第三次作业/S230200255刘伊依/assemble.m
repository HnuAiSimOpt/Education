function k_t=assemTriangle(k_t,k_ele,node1,node2,node3)
%assemTriangle This function assembles the element stiffness
% matrix k of the plane Triangle element with nodes
% i and j into the global stiffness matrix K.
% This function returns the global stiffness
% matrix K after the element stiffness matrix
% k is assembled.

d(1:2)=2*node1-1:2*node1;
d(3:4)=2*node2-1:2*node2;
d(5:6)=2*node3-1:2*node3;
for ii=1:6
    for jj=1:6
        k_t(d(ii),d(jj))=k_t(d(ii),d(jj))+k_ele(ii,jj);
    end
end
