%==========================================================================%
%                  HEAD                                                    %
%              TestPlatform                                                %
%          CAOYEHONG                                                       %
%              HUNAN UNIVERSITY.                                           %
%                                                          %
%==========================================================================%
function [gcoords,nodes,nnpe,ndpn,ElasticM,PoissonR,BC,load,dimension,element_type,SecType,Thickness,Density] = loadinp(dictionary)
format long
fid_inp= fopen(dictionary,'r');
temp_all = textscan( fid_inp ,'%s');
all=temp_all{1};
clear temp_all;
fclose(fid_inp);
clear fid_inp;
format long
keyword_node      = '*Node';
keyword_element   = '*Element,';
keyword_material  = '*Material,';
keyword_load      = '*Cload';
keyword_boundary  = '*Boundary';
keyword_sec       = 'Section:';
keyword_split1    = '*End';
keyword_split2    = 'Part';
location_split    = 0;
location_node     = 0;
location_sec      = 0;
location_element  = 0;
location_material = 0;
location_load     = [];
location_boundary = [];
bc_nodeset_temp   = [];
bctypedof         = [];
flag_load         = 0;
flag_boundary     = 0;
number_node       = 0;
number_element    = 0;
BC                = [];
bcset_smrk        = {};
bctype            = [];
gcoords           = zeros(1);
nodes             = zeros(1);
load              = zeros(1); 
Thickness         = 0;
SecType           = {};
count             = 1;
Density           = 0;
while location_split==0
    if strcmp(keyword_split1,all{count})&&strcmp(keyword_split2,all{count+1})
        location_split = count;
    end
    count=count+1;
end
for count = 1:location_split
    if strcmp(keyword_node,all{count})&&location_node==0
        location_node = count;
    end
    if strcmp(keyword_element,all{count})&&location_element==0
        location_element = count;
    end
    if strcmp(keyword_sec,all{count})
        location_sec = count;
    end
end
for count = location_split:length(all)
    if strcmp(keyword_material,all{count})
        location_material = count;
    end

    if strcmp(keyword_load,all{count})
        count2 = count;
        location_load(flag_load+1) = count;
        while ~strcmp('**',all{count2+4})
            flag_load = flag_load+1;
            location_load(flag_load+1) = count2+3;
            count2=count2+3;
        end
        flag_load=flag_load+1;
    end
    if strcmp(keyword_boundary,all{count})
        location_boundary(flag_boundary+1) = count;
        flag_boundary = flag_boundary+1;
    end
end
element_type = all{location_element+1};
temp=strsplit(all{location_element+1},'=');
element_type=temp{2};
fprintf('Element type:%s\n',element_type);
switch element_type
    case 'S4R' 
        nnpe=4;   
        ndpn=6;          
        dimension=3;    
    case 'S4R5'         
        nnpe=4;
        ndpn=5;
        dimension=3;
    case 'S3R'   
        nnpe=3;
        ndpn=6;
        dimension=3;
    case 'S4'    
        nnpe=4;
        ndpn=6;
        dimension=3;
    case 'S3'     
        nnpe=3;
        ndpn=6;
        dimension=3;
    case 'C3D4'  
        nnpe=4;
        ndpn=3;
        dimension=3;
    case 'C3D6'   
        nnpe=6;
        ndpn=3;
        dimension=3;
    case 'C3D8'   
        nnpe=8;
        ndpn=3;
        dimension=3;
    case 'C3D8R'  
        nnpe=8;
        ndpn=3;
        dimension=3;
    case 'C3D10M'    
        nnpe=10;
        ndpn=3;
        dimension=3;
    case 'CPS6'     
        nnpe=6;
        ndpn=2;
        dimension=2;
    case 'B31'
        nnpe=2;
        ndpn=2;
        dimension=3;
    case 'C3D8I'
        nnpe=8;
        ndpn=3;
        dimension=3;
    case 'CPS4'
        nnpe=4;
        ndpn=2;
        dimension=2;
    otherwise
        fprintf('Err:Element %s Not Supported\n',element_type);
end
count = location_node+1;
while count<location_element
    column = rem(count-location_node,dimension+1); 
    row = floor((count-location_node)/(dimension+1))+1;
    if column == 0
        row = row-1;
        column = dimension+1;
    end
    gcoords(row,column) = any2num(all{count},2);
    count = count+1;
end
count = location_element+2;
while str2num(all{count})
    column = rem(count-location_element-1,nnpe+1);
    row = floor((count-location_element-1)/(nnpe+1))+1;
    if column==0
        row=row-1;
        column=nnpe+1;
    end
    nodes(row,column)=any2num(all{count},1);
    count = count+1;
end
if strcmp('generate',all{count+2})
    num_end = any2num(all{count+4},2);
    num_start = any2num(all{count+3},2);
    num_gap = any2num(all{count+5},2);
    number_node = (num_end-num_start)/num_gap+1;
    count=count+6;
else
    while str2num(all{count+2})
        number_node = number_node+1;
        count = count+1;
    end
end
if length(gcoords(:,1))==number_node
    fprintf('Number of Nodes:%d\n',number_node);
else
    fprintf('Err:There was an error extracting nodes information!\n');
end
if strcmp('generate',all{count+2})
    num_end = any2num(all{count+4},2);
    num_start = any2num(all{count+3},2);
    num_gap = any2num(all{count+5},2);
    number_element = (num_end-num_start)/num_gap+1;
else
    while str2num(all{count+2})
        number_element = number_element+1;
        count = count+1;
    end
end
if length(nodes(:,1))==number_element
    fprintf('Number of Elements:%d\n',number_element);
else
    fprintf('Err:There was an error extracting elements information!\n');
end
if strcmp('*Density',all{location_material+2})
    Density = any2num(all{location_material+3},2);
    location_material=location_material+2;
end
if strcmp('*Elastic',all{location_material+2}) 
    ElasticM=any2num(all{location_material+3},2);
    PoissonR=any2num(all{location_material+4},2);
else
    fprintf('Err:Only one elastic material is supported!\n');
end
location_boundary_temp=[];
for count = 1:flag_boundary
    bcset_smrk{count} = strcat('nset=',all{location_boundary(count)+1});
    if str2num(all{location_boundary(count)+2}) 
        bctype{count} = 'dfbdof';
        location_boundary_temp(count) = location_boundary(count)+2;
    else
        bctype{count} = all{location_boundary(count)+2};
    end
end
for count = 1:length(all)
    for count2 = 1:length(bcset_smrk)
        if strcmp(bcset_smrk{count2},all{count}) 
            location_boundary(count2)=count;
        end
        count2=count2+1;
    end
    count=count+1;
end
flag=1;
for count=1:length(bctype)
    if strcmp('generate',all{location_boundary(count)+2})
        num_start = any2num(all{location_boundary(count)+3},2);
        num_end = any2num(all{location_boundary(count)+4},2);
        num_gap = any2num(all{location_boundary(count)+5},2);
        bc_nodeset_temp = num_start:num_gap:num_end;
    else
        count2 = location_boundary(count)+2;
        while str2num(all{count2})
            bc_nodeset_temp(flag) = any2num(all{count2},2);
            flag=flag+1;
            count2=count2+1;
        end
    end
    switch bctype{count}
        case 'XSYMM'
            bctypedof = [1,5,6];
        case 'YSYMM'
            bctypedof = [2,4,6];
        case 'ZSYMM'
            bctypedof = [3,4,5];
        case 'XASYMM'
            bctypedof = [2,3,4];
        case 'YASYMM'
            bctypedof = [2,4,6];
        case 'ZASYMM'
            bctypedof = [1,2,6];
        case 'PINNED'
            bctypedof = [1,2,3];
        case 'ENCASTRE'
            bctypedof = [1,2,3,4,5,6];
        case 'dfbdof' 
            while str2num(all{location_boundary_temp(count)}) 
                bctypedof = [bctypedof, str2num(all{location_boundary_temp(count)})];
                location_boundary_temp(count)=location_boundary_temp(count)+1;
                count=count+1;
            end
        otherwise
            fprintf('Err:This type of boundary constraint(%s) is not supported!\n',bctype{count});
    end
    count2=1;
    while count2<=length(bctypedof)
        if bctypedof(count2)>ndpn
            bctypedof(count2)=[];
        else
            count2=count2+1;
        end
    end
    for count2 = 1:length(bc_nodeset_temp)
        BC=[BC (bc_nodeset_temp(count2)-1)*ndpn+bctypedof];
        count2=count2+1;
    end
end
fprintf('The total number constrained DOFs is :%d\n',length(BC));
while ~contains(all{location_sec},'*')
    location_sec=location_sec+1;
end
temp=strsplit(all{location_sec},'*');
SecType=temp{2};
count2=location_sec;
if strcmp(all{count2},'*Shell')
   while count2<location_sec+100
       while str2num(all{count2})
           Thickness = any2num(all{count2},2);
           count2 = location_sec+100;
           break;
       end
       count2=count2+1;
   end
end
keyword_load={}; 
for count = 1:length(location_load)
    keyword_load{count}=strcat('nset=',all{location_load(count)+1});
end
location_load_temp = location_load+1;
flag_load = 1;
for count = location_split:length(all)
    for count2=1:length(keyword_load)
        if strcmp(keyword_load{count2},all{count})
            location_load(flag_load) = count; 
            flag_load = flag_load+1;
        end
    end
end
flag=1;
for count = 1:length(location_load)
    location_temp=location_load(count)+2;
    while any2num(all{location_temp},2)
        load_dof = any2num(all{location_load_temp(count)+1},1);
        load_node = any2num(all{location_temp},1);
        load(flag,1)= (load_node-1)*ndpn+load_dof;
        load(flag,2) = any2num(all{location_load_temp(count)+2},2);
        location_temp=location_temp+1;
        flag=flag+1;
    end
end
end