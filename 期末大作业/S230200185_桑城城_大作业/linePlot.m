function linePlot(xyzi,vari,min_vari,max_vari)
%%%%% 绘图
x = xyzi(:,1); 
y = xyzi(:,2);
z = xyzi(:,3);

surface([x,x], [y,y], [z,z], [vari,vari],...
	'FaceColor', 'no',...
	'EdgeColor', 'interp',...
	'LineWidth', 4);
grid on;
colormap jet;
caxis([min_vari max_vari]);
h = colorbar;
zstep=(max_vari-min_vari)/10;
set(h,'ytick',(max_vari:zstep:max_vari)');
end

