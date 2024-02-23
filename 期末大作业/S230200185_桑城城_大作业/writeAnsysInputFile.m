format long;
formatSpec = '%.15f';

fileID = fopen('input_load_ANSYS.txt','w');
fprintf(fileID,'\n');
fprintf(fileID,'\n');

fprintf(fileID,'/SOLU\n');

fprintf(fileID,'!!!!!!**Finite Element Analysis with MATLAB and ANSYS**!!!!!\n');
fprintf(fileID,'!!!!!!*******By: CONG-NGUYEN******!!!!!\n');
fprintf(fileID,'!!!!!!*******3D Building Frames******!!!!!\n');
fprintf(fileID,'\n');


free_loadP = setdiff((1:nP)', loadP);
loadP = setdiff((1:nP)', free_loadP);

for i = 1:length(loadP)
   nodei = loadP(i,1);
   vari = ['F,',num2str(nodei),',FX,',num2str(force(nodei,1),formatSpec),'\n'];
   fprintf(fileID,vari);

end

fprintf(fileID,'ALLSEL\n');
fprintf(fileID,'\n');

fprintf(fileID,'SOLVE\n');
fprintf(fileID,'FINISH \n');
fclose(fileID);

