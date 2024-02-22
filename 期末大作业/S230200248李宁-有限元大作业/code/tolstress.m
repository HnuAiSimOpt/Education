function stress_node=tolstress(nel,nnode,nodes,stress)

neigh_node = cell(nnode,1);
neigh_node_ind = cell(nnode,1);
indneigh=zeros(1,nnode);
for i=1:nel
    for j=1:4
        indneigh(nodes(i,j))=indneigh(nodes(i,j))+1;
        neigh_node{nodes(i,j)}(indneigh(nodes(i,j)))=i;
        neigh_node_ind{nodes(i,j)}(indneigh(nodes(i,j)))=j;
    end
end

stress_node=zeros(3,nnode);	
for inode=1:nnode
    numel= indneigh(inode);
    for i=1:numel
        ind_nel= neigh_node{inode}(i);
        ind_nod=neigh_node_ind{inode}(i);
        for j=1:3
            stress_node(j,inode)=stress_node(j,inode)+stress(ind_nel,ind_nod,j);
        end
    end
    stress_node(:,inode)=stress_node(:,inode)/numel;
end
