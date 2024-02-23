function scatterplot2(xyz,phi,phimin,phimax,sizescatter)
tol = 1e-4;
c1 = zeros(3,1);c2 = ones(3,1);c33 = 0:1/2:1;c3 = c33';c44 = 1:-1/2:0;
c55 = 1/2:-1/2:0;c4 = c44';cc1 = [c1,c3,c4];cc2 = [c3,c2,c1];
cc3 = [ones(2,1),c55',zeros(2,1)];
cmap = [cc1;cc2(2:end,:);cc3];


colormap(cmap);
ncolor = size(cmap,1);
ztepz = (phimax-phimin)/ncolor;
colors = zeros(1,length(phi));
colors(phi<= phimin+tol)= phimin;
colors(phi> phimin+tol & phi <= phimin+ztepz)= phimin + 1*ztepz;
colors(phi>phimin+ztepz & phi<=phimin+2*ztepz)= phimin+2*ztepz;
colors(phi>phimin+2*ztepz & phi<=phimin+3*ztepz)= phimin+3*ztepz;
colors(phi>phimin+3*ztepz & phi<=phimin+4*ztepz)= phimin+4*ztepz;
colors(phi>phimin+4*ztepz & phi<=phimin+5*ztepz)= phimin+5*ztepz;
colors(phi>phimin+5*ztepz & phi<=phimin+6*ztepz)= phimin+6*ztepz;
colors(phi>phimin+6*ztepz)= phimax;

scatter3(xyz(:,1),xyz(:,2),xyz(:,3),sizescatter,colors,'filled');
caxis([phimin phimax]);
h = colorbar;
tickrange = (phimin:ztepz:phimax);
TL=arrayfun(@(x) sprintf('%.3f',x),tickrange,'un',0);
set(h,'ytick',tickrange);
set(h,'TickLabels',TL)

end