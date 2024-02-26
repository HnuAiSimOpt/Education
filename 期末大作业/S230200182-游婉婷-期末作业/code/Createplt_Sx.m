function []=Createplt_Sx(node,element_solid,sx)
eles_num=size(element_solid,1);
nodes_num=size(node,1);
for i=1:1860
    Sx(element_solid(i,2))=sx(i);
    Sx(element_solid(i,3))=sx(i);
    Sx(element_solid(i,4))=sx(i);
    Sx(element_solid(i,5))=sx(i);
    Sx(element_solid(i,6))=sx(i);
    Sx(element_solid(i,7))=sx(i);
    Sx(element_solid(i,8))=sx(i);
    Sx(element_solid(i,9))=sx(i);
end
fid=fopen('result_Sx.plt','w');
fprintf(fid,'TITLE="test case governed by poisson equation"\n');
fprintf(fid,'VARIABLES="x""y""z""sigax"\n');
fprintf(fid,'ZONE N=%8d,E=%8d,ET=BRICK,F=FEPOINT\n',nodes_num,eles_num);
for i = 1:nodes_num
    fprintf(fid,'%16.6e%16.6e%16.6e%16.6e\n',node(i,2),node(i,3),node(i,4),Sx(1,i));
end
for i=1:eles_num
    fprintf(fid,'%8d%8d%8d%8d%8d%8d%8d%8d\n',element_solid(i,2),element_solid(i,3),element_solid(i,4),element_solid(i,5),element_solid(i,6),element_solid(i,7),element_solid(i,8),element_solid(i,9));
end
fclose(fid);