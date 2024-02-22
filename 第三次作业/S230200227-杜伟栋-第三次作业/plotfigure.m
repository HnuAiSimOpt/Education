<<<<<<< HEAD
function plotfigure(element,node,linecolor)
% linecolor='k--'
numelem=size(element,1);
for iel=1:numelem
    nod=element(iel,:);
    plot(node(nod,1),node(nod,2),linecolor)
    hold on;
end
=======
function plotfigure(element,node,linecolor)
% linecolor='k--'
numelem=size(element,1);
for iel=1:numelem
    nod=element(iel,:);
    plot(node(nod,1),node(nod,2),linecolor)
    hold on;
end
>>>>>>> 0c2441c35348cc2cb0d4bfa4f93dc736c7d002fd
