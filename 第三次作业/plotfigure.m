function plotfigure(element,node,linecolor)
% linecolor='k--'
numelem=size(element,1);
for iel=1:numelem
    nod=element(iel,:);
    plot(node(nod,1),node(nod,2),linecolor)
    hold on;
end
