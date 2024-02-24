function scatterplot1(xyz,vari,varimin,varimax,sizescatter)
c1 = zeros(6,1);c2 = ones(6,1);c33 = 0:0.2:1;c3 = c33';c44 = 1:-0.2:0;
c55 = 0.8:-0.2:0;c4 = c44';cc1 = [c1,c3,c4];cc2 = [c3,c2,c1];
cc3 = [ones(5,1),c55',zeros(5,1)];cmap = [cc1;cc2;cc3];

colormap(cmap);
ncolor = size(cmap,1);
ztepzscaled = (varimax-varimin)/ncolor;
ztepz = (varimax-varimin)/ncolor;
colors = zeros(1,length(vari));
colors(vari==varimin)=varimin;
colors(vari<varimin+ztepzscaled)=varimin + 1*ztepz;
colors(vari>=varimin+ztepzscaled & vari<varimin+2*ztepzscaled)=varimin+2*ztepz;
colors(vari>=varimin+2*ztepzscaled & vari<varimin+3*ztepzscaled)=varimin+3*ztepz;
colors(vari>=varimin+3*ztepzscaled & vari<varimin+4*ztepzscaled)=varimin+4*ztepz;
colors(vari>=varimin+4*ztepzscaled & vari<varimin+5*ztepzscaled)=varimin+5*ztepz;
colors(vari>=varimin+5*ztepzscaled & vari<varimin+6*ztepzscaled)=varimin+6*ztepz;
colors(vari>=varimin+6*ztepzscaled & vari<varimin+7*ztepzscaled)=varimin+7*ztepz;
colors(vari>=varimin+7*ztepzscaled & vari<varimin+8*ztepzscaled)=varimin+8*ztepz;
colors(vari>=varimin+8*ztepzscaled & vari<varimin+9*ztepzscaled)=varimin+9*ztepz;
colors(vari>=varimin+9*ztepzscaled & vari<varimin+10*ztepzscaled)=varimin+10*ztepz;
colors(vari>=varimin+10*ztepzscaled & vari<varimin+11*ztepzscaled)=varimin+11*ztepz;
colors(vari>=varimin+11*ztepzscaled & vari<varimin+12*ztepzscaled)=varimin+12*ztepz;
colors(vari>=varimin+12*ztepzscaled & vari<varimin+13*ztepzscaled)=varimin+13*ztepz;
colors(vari>=varimin+13*ztepzscaled & vari<varimin+14*ztepzscaled)=varimin+14*ztepz;
colors(vari>=varimin+14*ztepzscaled & vari<varimin+15*ztepzscaled)=varimin+15*ztepz;
colors(vari>=varimin+15*ztepzscaled )=varimax;
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),sizescatter,colors,'filled');
caxis([varimin varimax]);
h = colorbar;
set(h,'ytick',varimin:ztepz:varimax);
end