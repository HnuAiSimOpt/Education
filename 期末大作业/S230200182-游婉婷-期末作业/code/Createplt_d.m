function []=Createplt_d(node,element_solid,d)
eles_num=size(element_solid,1);
nodes_num=size(node,1);
fid=fopen('result_d.plt','w');   % 以只写模式打开result_d.plt文件
fprintf(fid,'TITLE="test case governed by poisson equation"\n');
fprintf(fid,'VARIABLES="x""y""z""U""V""W"\n');
fprintf(fid,'ZONE N=%8d,E=%8d,ET=BRICK,F=FEPOINT\n',nodes_num,eles_num);
for i = 1:nodes_num
fprintf(fid,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',node(i,2),node(i,3),node(i,4),d(3*i-2),d(3*i-1),d(3*i));   % 输出节点坐标及对应方向位移
end
for i=1:eles_num
fprintf(fid,'%8d%8d%8d%8d%8d%8d%8d%8d\n',element_solid(i,2),element_solid(i,3),element_solid(i,4),element_solid(i,5),element_solid(i,6),element_solid(i,7),element_solid(i,8),element_solid(i,9));   % 输出单元对应8个节点编号
end
fclose(fid);