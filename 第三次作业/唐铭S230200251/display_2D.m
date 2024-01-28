function display_2D(elemnode,nodex,nodey,nodexlocold,nodeylocold)

face = [1 2 3];
    set(gcf,'Name','ISO display','NumberTitle','off');
    for i = 1:size(elemnode)
                    vert1=[nodex(elemnode(i,:)) nodey(elemnode(i,:))];
                    patch('Faces',face,'Vertices',vert1,'FaceColor','none','LineStyle','-','EdgeColor','r');
                    vert2=[nodexlocold(elemnode(i,:)) nodeylocold(elemnode(i,:))];
                    patch('Faces',face,'Vertices',vert2,'FaceColor','none','LineStyle','--','EdgeColor','b');
                    hold on;
    end
    axis equal; 
    axis tight; 
    axis off; 
    box on;
    pause(1e-6);
end