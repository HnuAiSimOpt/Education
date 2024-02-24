%==========================================================================%
%                  HEAD                                                    %
%              TestPlatform_explicit                                       %
%          CAOYEHONG                                                       %
%              HUNAN UNIVERSITY.                                           %
%                                                                          %
%==========================================================================%

% START
clear;
close all;
clc;
format long
%==========================================================================

% LOAD DATA
addpath('Results_Model');
fprintf('可计算单元类型:\nCPS4-线性平面四边形单元\nC3D4-线性四面体单元\n');
dictionary = input('请输入单元名：\n如CPS4\n','s');
dictionary = strcat(dictionary,'.inp');
[gcoords,nodes,nnpe,ndpn,ElasticM,PoissonR,BC,load,dimension,Element_Type,SecType,Thickness,Density] = loadinp(dictionary);
fprintf('Load Data Finished\n');
%==========================================================================

% PAPERING WORK
     sT = 0;                                %开始时间
     dT = 0.05;    % 时间步长
     eT = 20;                                %结束时间
  stepn = floor((eT-sT)/dT);
maxloop = 12;
    Err = 100;
      C = 1;  % 阻尼
    nes = length(nodes(:,1));   
    nns = length(gcoords(:,1));   
    nds = nns * ndpn;              
   ndpe = ndpn * nnpe;                
     kk = zeros(nds,nds);
    fin = zeros(nds,stepn);            
    fex = zeros(nds,1);                
     ff = zeros(nds,1); 
   disp = zeros(nds,1);
    vel = zeros(nds,1);
    acc = zeros(nds,1);
  index = sparse(ndpe,1);            
%==========================================================================


% CONSTITUTIVE MODELS
switch Element_Type
    case 'C3D4'
        Cm = fematiso(4,ElasticM,PoissonR);
    case 'CPS4'
        Cm = fematiso(1,ElasticM,PoissonR);
    otherwise
        fprintf('NOT SUPPORTED YET\n');
        return;
end
%==========================================================================
for count_load=1:length(load(:,1))
    fex(load(count_load,1))=load(count_load,2);
end

% COMPUTE MASS MATRIX
switch Element_Type
    case 'C3D4'
        Mass = calmass_C3D4(gcoords,nodes,Density,Thickness,nnpe,nds);
        kk=getK_C3D4(gcoords,nodes,nes,Cm,nnpe,ndpn,kk);
    case 'CPS4'
        Mass = calmass_CPS4(gcoords,nodes,Density,1,nnpe,nds);
        kk=getK_CPS4(gcoords,nodes,nes,Cm,nnpe,ndpn,kk);
    otherwise
        fprintf('NOT SUPPORTED YET\n');
        return;
end

invm=inv(Mass);
countloop=1;
flagloop=0;
while countloop<=maxloop && flagloop==0
    countstep=1;
    flagstep=1;
    while countstep<=stepn && flagstep==1
        fin(:,countstep)=kk*disp(:,countstep);
        ff(:,countstep)=fex(:,1)-fin(:,countstep)-C*vel(:,countstep);
        acc(:,countstep)=invm*ff(:,countstep);
        for count1=1:length(BC)
            acc(BC(count1),:)=0;
        end
        vel(:,countstep+1)=vel(:,countstep)+acc(:,countstep)*dT;
        disp(:,countstep+1)=disp(:,countstep)+vel(:,countstep+1)*dT;
        if max(disp(:,countstep+1))>Err
            fprintf('dT = %f does not converge\n',dT);
            flagstep=0;
        end
        countstep=countstep+1;
    end
    if flagstep==1
        flagloop=1;
    else
        dT=dT/2;
        fprintf('The time step is halved\ndT = %f\n\n',dT);
        stepn = floor((eT-sT)/dT);
    end
    countloop=countloop+1;
    if countloop>maxloop
        fprintf('Err:The results do not converge\n');
        return;
    end
end
%==========================================================================

result='result';
cd('Results_Model\ExplicitResults');% get into result file.
fprintf('OUTPUTing...\n');
for countout=1:stepn
    name=[result,num2str(countout)];
    name=strcat(name,'.plt');
    fid = fopen(name,'w');
    switch Element_Type
        case 'CPS4'
            fprintf(fid,'TITLE="CPS4 EXPLICIT"\n');
            fprintf(fid,'VARIABLES="X""Y""DispX""DispY"\n');
            fprintf(fid,'ZONE T="flow-field",N=%8d,E=%8d,ET=QUADRILATERAL,F=FEPOINT\n',nns,nes);
            for count = 1:nns
                fprintf(fid,'%16.6e %16.6e %16.6e %16.6e\n',gcoords(count,2),gcoords(count,3),disp(2*count-1,countout),disp(2*count,countout));
            end
            for count=1:nes
                fprintf(fid,'%8d%8d%8d%8d\n',nodes(count,2),nodes(count,3),nodes(count,4),nodes(count,5));
            end
        case 'C3D4'
            fprintf(fid,'TITLE="C3D4 EXPLICIT"\n');
            fprintf(fid,'VARIABLES="X""Y""Z""DispX""DispY""DispZ"\n');
            fprintf(fid,'ZONE T="flow-field",N=%8d,E=%8d,ET=TETRAHEDRON,F=FEPOINT\n',nns,nes);
            for count = 1:nns
                fprintf(fid,'%16.6e %16.6e %16.6e %16.6e %16.6e %16.6e\n',gcoords(count,2),gcoords(count,3),gcoords(count,4),disp(3*count-2,countout),disp(3*count-1,countout),disp(3*count,countout));
            end
            for count=1:nes
                fprintf(fid,'%8d%8d%8d%8d\n',nodes(count,2),nodes(count,3),nodes(count,4),nodes(count,5));
            end
        otherwise
            return;
    end
    fclose(fid);
end
fprintf('total number of files:%d\n',stepn);
cd('..\..');
%=========================================================================%
%                                END CODE                                 %
%=========================================================================%
