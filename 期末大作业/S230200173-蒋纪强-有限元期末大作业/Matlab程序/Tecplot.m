function  Tecplot(IEN,gcoord,stress_n)
% 将结果导出到Mise.plt程序文件
nel = size(IEN,1);
node= size(gcoord,1);

nodedata=fopen('stress.plt','w');  
fprintf(nodedata,'TITLE="data"\n');
fprintf(nodedata,'VARIABLES=,"X", "Y","Z" ,"sig1","sig2","sig3","sig4","sig5","sig6","sig7"\n');
fprintf(nodedata,'ZONE T="%d  "  ,  N=%d, E=%d, ET=BRICK, F=FEPOINT\n',1,node,nel);

sigm1=stress_n(:,1);
sigm2=stress_n(:,2);
sigm3=stress_n(:,3);
sigm4=stress_n(:,4);
sigm5=stress_n(:,5);
sigm6=stress_n(:,6);
sigm7=stress_n(:,7);

for i=1:node
    fprintf(nodedata,'%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f,%10.10f\n',...
      gcoord(i,1),gcoord(i,2),gcoord(i,3),sigm1(i),sigm2(i),sigm3(i),sigm4(i),sigm5(i),sigm6(i),sigm7(i));
end  

for i=1:nel
    for j=1:size(IEN,2)
        fprintf(nodedata,'%d       ',IEN(i,j));  
    end
    fprintf(nodedata,'\n');
end
end
