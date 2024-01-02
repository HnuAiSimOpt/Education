%--------------------------------------------------------
%purpose:
%find the connected node and coord values of each element
%
%Synopsis:
%[nd,x,y]=cal_nd_x_y(iel,nodes,gcoord)
%
%variable description:
%nd_x_y - the connected node and coord values of each element
%iel - iel th element 
%nodes - data for nodal connectivity for each element
%gcoord - the coord data of all element
%--------------------------------------------------------

function [nd,x,y]=cal_nd_x_y(iel,nodes,gcoord)

nsel=length(nodes(iel,:));

for ie=1:nsel
    
    nd(ie)=nodes(iel,ie);         %ie th connected node for (iel)-th element
    
    x(ie)=gcoord(nd(ie),1); y(ie)=gcoord(nd(ie),2);     %coord values of ie th node

end