function structural_stiffness_matrix=TetrahedronAssemble(element_stiffness_matrix,element_node,node_coordinate)
[element_number,n_ele]=size(element_node);
node_number=size(node_coordinate,1);
structural_stiffness_matrix=zeros(node_number*3,node_number*3);
for i=1:element_number
    for m=1:n_ele
        for n=1:n_ele
            i_m=element_node(i,m);
            i_n=element_node(i,n);
            structural_stiffness_matrix(3*i_m-2:3*i_m,3*i_n-2:3*i_n)=structural_stiffness_matrix(3*i_m-2:3*i_m,3*i_n-2:3*i_n)+element_stiffness_matrix(3*m-2:3*m,3*n-2:3*n);
        end
    end
end