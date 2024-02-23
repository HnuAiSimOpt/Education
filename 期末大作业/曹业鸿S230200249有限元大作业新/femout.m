function []=femout(Disp,gcoords,nodes,nns,ndpn,Element_Type)
[x,y,z,M,tx,ty,tz,Mt] = DataPaD(Disp,nns,ndpn);
nes=length(nodes(:,1));
name=['result_',Element_Type];
name=[name,'.plt'];
switch Element_Type
    case 'CPS6'
        flag=1;
    case 'CPS4'
        flag=2;
    case 'C3D4'
        flag=3;
    case 'S4R'
        flag=5;
    otherwise
        fprintf('FEMout Wrong\n');
end
cd('Results_Model');
fid = fopen(name,'w');
fprintf(fid,'TITLE="%s"\n',Element_Type);
switch flag
    
    case 1 % CPS6--6节点 2自由度 平面
        fprintf(fid,'VARIABLES="X""Y""U""V""M"\n');
        maxnode=max([max(nodes(:,2)),max(nodes(:,3)),max(nodes(:,4))]);
        fprintf(fid,'ZONE T="flow-field",N=%8d,E=%8d,ET=TRIANGLE,F=FEPOINT\n',maxnode,nes);
        for count = 1:maxnode
            fprintf(fid,'%16.6e %16.6e %16.6e %16.6e %16.6e\n',gcoords(count,2),gcoords(count,3),x(count),y(count),M(count));
        end
        for count=1:nes
            fprintf(fid,'%8d %8d %8d\n',nodes(count,2),nodes(count,3),nodes(count,4));
        end
        
    case 2 % CPS4--4节点 2自由度 平面
        fprintf(fid,'VARIABLES="X""Y""DispX""DispY""Mdisp"\n');
        fprintf(fid,'ZONE T="flow-field",N=%8d,E=%8d,ET=QUADRILATERAL,F=FEPOINT\n',nns,nes);
        for count = 1:nns
            fprintf(fid,'%16.6e %16.6e %16.6e %16.6e %16.6e\n',gcoords(count,2),gcoords(count,3),x(count),y(count),M(count));
        end
        for count=1:nes
            fprintf(fid,'%8d%8d%8d%8d\n',nodes(count,2),nodes(count,3),nodes(count,4),nodes(count,5));
        end
        
    case 3 % C3D4--4节点 3自由度 空间
        fprintf(fid,'VARIABLES="X""Y""Z""U""V""W""M"\n');
        fprintf(fid,'ZONE T="flow-field",N=%8d,E=%8d,ET=TETRAHEDRON,F=FEPOINT\n',nns,nes);
        for count = 1:nns
            fprintf(fid,'%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e%16.6e\n',gcoords(count,2),gcoords(count,3),gcoords(count,4),x(count),y(count),z(count),M(count));
        end
        for i=1:nes
            fprintf(fid,'%8d %8d %8d %8d\n',nodes(i,2),nodes(i,3),nodes(i,4),nodes(i,5));
        end
        
    case 5 % S4R --4节点 6自由度 空间
        fprintf(fid,'VARIABLES="X""Y""Z""DispX""DispY""DispZ""Mdisp""thetaX""thetaY""thetaZ""Mtheta"\n');
        fprintf(fid,'ZONE T="flow-field",N=%8d,E=%8d,ET=QUADRILATERAL,F=FEPOINT\n',nns,nes);
        for count = 1:nns
            fprintf(fid,'%16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e %16.6e\n',gcoords(count,2),gcoords(count,3),gcoords(count,4),x(count),y(count),z(count),M(count),tx(count),ty(count),tz(count),Mt(count));
        end
        for count=1:nes
            fprintf(fid,'%8d%8d%8d%8d\n',nodes(count,2),nodes(count,3),nodes(count,4),nodes(count,5));
        end
        
    otherwise
        fprintf('FEMout Wrong\n');
end

fclose(fid);
cd ..;
return