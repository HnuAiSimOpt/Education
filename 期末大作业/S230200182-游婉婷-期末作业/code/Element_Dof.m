function [index]=Element_Dof(nod,nel,node_dof)
   k=0;
   for i=1:nel     
     temp=(nod(i)-1)*node_dof;
       for j=1:node_dof   
         k=k+1;
         index(k)=temp+j;
       end
   end

