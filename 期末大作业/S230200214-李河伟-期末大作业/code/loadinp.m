
function [gcoords,nodes,nnpe,ndpn,ElasticM,PoissonR,BC,load,dimension,element_type,SecType,Thickness,Density] = loadinp(dictionary)
format long
%% Extract the info from INP exported from abaqus.
%      gcoords -- �ڵ�������Ϣ [�ڵ���,�ڵ�����1,�ڵ�����2,�ڵ�����3...;]
%     ElasticM -- ����ģ��
%     PoissonR -- ���ɱ�
%        nodes -- ��Ԫ������Ϣ [��Ԫ���,�ڵ�1,�ڵ�2,�ڵ�3...;]
%           BC -- �߽����������ɶȣ�   Ĭ�ϱ߽�����Լ��Ϊ0
%         load -- �غ� [�����غɵ����ɶ�,�غɴ�С;]
%         nnpe -- ��Ԫ�ڵ������뵥Ԫ�����й�
%         ndpn -- �ڵ����ɶ���
% all���Կո�Ϊ�ָ��Ķ�ȡ�ַ��������е���Ϣ��ȡΪһ��ϸ���������all�У�ͨ��all{?}����

%% �ļ�׼��---��ȡ
fid_inp= fopen(dictionary,'r');
temp_all = textscan( fid_inp ,'%s');
all=temp_all{1};    % ��Ƕ��
clear temp_all;     % �ͷ��ڴ�
fclose(fid_inp);    % �ر��ļ�
clear fid_inp;
format long

%% �ؼ��ּ��ڴ����
keyword_node      = '*Node';         % node��ʼ�ı�ǡ����ο�inp�ļ��淶
keyword_element   = '*Element,';     % ��Ԫ��ʼ
keyword_material  = '*Material,';    % ���ϵĹؼ���
keyword_load      = '*Cload';        % �����غɹؼ���
keyword_boundary  = '*Boundary';     % Լ���Ĺؼ���
keyword_sec       = 'Section:';
keyword_split1    = '*End';          % ��ǰ��part������Լ��֮��ķֿ�
keyword_split2    = 'Part';
location_split    = 0;
location_node     = 0;
location_sec      = 0;
location_element  = 0;
location_material = 0;
location_load     = [];
location_boundary = [];
bc_nodeset_temp   = [];              % Լ���㼯
bctypedof         = [];              % Լ�����Ͷ�ӦԼ�������ɶ�
flag_load         = 0;               % �����ж��ּ����غ�
flag_boundary     = 0;               % �����ж��Լ������
number_node       = 0;               % inp�ļ��м�¼�Ľڵ���
number_element    = 0;               % inp�ļ��м�¼�ĵ�Ԫ�������ڶ�ȡ���ĵ�Ԫ�����к˶�
BC                = [];              % ��Լ�������ɶ�
bcset_smrk        = {};              %  
bctype            = [];              % Լ�����ͣ�֧��abaqusԼ����ȫ������
gcoords           = zeros(1);
nodes             = zeros(1);
load              = zeros(1);        % [�غ��������ɶ�1,�غɴ�С1;
                                     %  �غ��������ɶ�2,�غɴ�С2;]
Thickness         = 0;               % SecType==solid --> Thickness=0
SecType           = {};
count             = 1;
Density           = 0;

%% ���ҹؼ��ֵ�λ��
while location_split==0
    if strcmp(keyword_split1,all{count})&&strcmp(keyword_split2,all{count+1}) % ���ļ��ֶ�
        location_split = count;
    end
    count=count+1;
end
for count = 1:location_split
    if strcmp(keyword_node,all{count})&&location_node==0  %����node��ʼ��λ��
        location_node = count;   %�ڵ����꿪ʼ��λ��
    end
    if strcmp(keyword_element,all{count})&&location_element==0
        location_element = count;     %��Ԫ������Ϣ��ʼ�ĵط�
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
        location_load(flag_load+1) = count;   % ��������غɻ���һ���ڵ�ͬʱʩ�Ӷ�������غ�
        while ~strcmp('**',all{count2+4})
            flag_load = flag_load+1;
            location_load(flag_load+1) = count2+3;
            count2=count2+3;
        end
        flag_load=flag_load+1;
    end
    if strcmp(keyword_boundary,all{count})
        location_boundary(flag_boundary+1) = count;   % ���Լ��
        flag_boundary = flag_boundary+1;
    end
end

%% ��Ԫ����---���ݵ�Ԫ����ȷ����Ԫ�ڵ��������ɶȣ���һ�����ܶ�ȡ
element_type = all{location_element+1};
temp=strsplit(all{location_element+1},'=');
element_type=temp{2};
fprintf('Element type:%s\n',element_type);
switch element_type
    case 'S4R'             % 4�ڵ��ı����������ֿǵ�Ԫ
        nnpe=4;            % ��Ԫ�ڵ���
        ndpn=6;            % �ڵ����ɶ���
        dimension=3;       % ��ά��Ԫ -- ȷ��ÿ���ڵ����������xy or xyz
    case 'S4R5'            % 4�ڵ��ı����������ֱ��ǵ�Ԫ��ɳ©����
        nnpe=4;
        ndpn=5;
        dimension=3;
    case 'type=S3R'        % 3�ڵ��������������ֿǵ�Ԫ
        nnpe=3;
        ndpn=6;
        dimension=3;
    case 'S4'              % 4�ڵ�ͨ�ÿǵ�Ԫ ����ĤӦ��
        nnpe=4;
        ndpn=6;
        dimension=3;
    case 'S3'              % 3�ڵ������οǵ�Ԫ ����ĤӦ��
        nnpe=3;
        ndpn=6;
        dimension=3;
    case 'C3D4'            % 4�ڵ����������Ե�Ԫ
        nnpe=4;
        ndpn=3;
        dimension=3;
    case 'C3D6'            % 6�ڵ�����������������Ԫ
        nnpe=6;
        ndpn=3;
        dimension=3;
    case 'C3D8'            % 8�ڵ�����ʵ�嵥Ԫ
        nnpe=8;
        ndpn=3;
        dimension=3;
    case 'C3D8R'           % 8�ڵ�����ʵ�嵥Ԫ �������� ɳ©����
        nnpe=8;
        ndpn=3;
        dimension=3;
    case 'C3D10M'          % 10�ڵ�������������ε�Ԫ ɳ©����
        nnpe=10;
        ndpn=3;
        dimension=3;
    case 'CPS6'            % 6�ڵ������ζ���ƽ��Ӧ����Ԫ
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

%% ��ȡ�ڵ��������
count = location_node+1;    % �ڵ��������ʼλ���ǹؼ��ֵ���һλ
while count<location_element  % �ڵ�����֮��������ǵ�Ԫ�Ĺؼ���
    column = rem(count-location_node,dimension+1);  % �У���ǰλ�ö�ά������һȡ�࣬�ڵ���-�ڵ�����1-�ڵ�����2-�ڵ�����3
    row = floor((count-location_node)/(dimension+1))+1; % �� ����ȡ����+1
    if column == 0  % һ�е����һλ������Ϊ0
        row = row-1;
        column = dimension+1;
    end
    gcoords(row,column) = any2num(all{count},2);   %��ȡС��
    count = count+1;
end

%% ��ȡ��Ԫ�ڵ�������Ϣ
count = location_element+2;% ��Ԫ�ؼ��� - ��Ԫ���� - ��Ԫ��Ϣ��ʼ   �ؼ��ֵ���Ԫ��Ϣ��ʼҪ+2
while str2num(all{count})  % ��Ԫ����Լ����ӵĽڵ���ȫΪ���������ڵ�Ԫ��Ϣ������Ϊ�ַ�����str2num֮��Ϊfalse���ڵ��������Ϊ0
    column = rem(count-location_element-1,nnpe+1);
    row = floor((count-location_element-1)/(nnpe+1))+1;
    if column==0
        row=row-1;
        column=nnpe+1;
    end
    nodes(row,column)=any2num(all{count},1);    %��ȡ����
    count = count+1;
end

%% ��֤�ڵ���Ŀ�͵�Ԫ��Ŀ��ȡ��ȷ
if strcmp('generate',all{count+2})  %�����countָ��Ԫ������ĵ�һ��λ��--��Ԫ������*Nset
    num_end = any2num(all{count+4},2); %�����count��������ȡ��Ԫ��Ϣ��count��������
    num_start = any2num(all{count+3},2);
    num_gap = any2num(all{count+5},2);
    number_node = (num_end-num_start)/num_gap+1;
    count=count+6;   % ����count+6��ָ���ע��Ԫ����λ�ã���*Elset��count���������ȡ��Ԫ��������
else
    while str2num(all{count+2})  %���ǵ�Ԫ����˳�����ɵ����������countһֱ���ƶ������ͣ��*Elset
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

%% ���Ҳ��ϲ����������Բ��ϡ���ֻ��һ�ֲ���
if strcmp('*Density',all{location_material+2})  % �����ڵ��Բ���
    Density = any2num(all{location_material+3},2);
    location_material=location_material+2;
end
if strcmp('*Elastic',all{location_material+2})  % �����ڵ��Բ���
    ElasticM=any2num(all{location_material+3},2);   % ��ȡС��
    PoissonR=any2num(all{location_material+4},2);   % ��ȡС��
else
    fprintf('Err:Only one elastic material is supported!\n');
end


%% ��ȡ�߽�����   ֧��ABAQUS���ֶ���Լ���ķ�ʽ����ȡ�߽�����ģ��������
% ���ҵ��߽�������ʼ��λ�ã��ҳ�����ĵ㼯��Լ�����ͣ�����Լ�����ͷ�Ϊ����Լ���Ĺؼ��ֶ����ֱ�Ӹ������ɶ�������
location_boundary_temp=[];
for count = 1:flag_boundary  %����countָʾԼ������Ŀ  �������ҵ�Լ���ĵ㼯����bcset_smrk�����Ӧ��Լ������bctype
    bcset_smrk{count} = strcat('nset=',all{location_boundary(count)+1});  % �߽�Լ���㼯
    if str2num(all{location_boundary(count)+2})   % ������λ�������ֵĻ���˵���ǵڶ��ֶ���Լ���ķ�ʽ
        bctype{count} = 'dfbdof';   % �ڶ��ֶ���Լ�������ͣ�ֱ�Ӹ���Լ�������ɶȣ��������Լ�Լ����һ���ؼ��֣��൱����һ��Լ������
        location_boundary_temp(count) = location_boundary(count)+2; %��¼������ֱ�Ӹ�����Լ�������ɶȵ�λ�ã�����location_boundary�����Լ���ĵ㼯��λ����
    else
        bctype{count} = all{location_boundary(count)+2};% Լ������  bctype������bcset_smrk,location_boundary�ȶ�Ӧ
    end
end

for count = 1:length(all)
    for count2 = 1:length(bcset_smrk)
        if strcmp(bcset_smrk{count2},all{count})    % ����Լ���㼯��λ��
            location_boundary(count2)=count;       % ����location_boundary�����Լ���㼯��λ��
        end
        count2=count2+1;
    end
    count=count+1;
end
flag=1;
for count=1:length(bctype)   %length(bctype)==length(location_boundary)  ����Լ�������ɶ�
    if strcmp('generate',all{location_boundary(count)+2})   % Լ���ĵ㼯��������
        num_start = any2num(all{location_boundary(count)+3},2);
        num_end = any2num(all{location_boundary(count)+4},2);
        num_gap = any2num(all{location_boundary(count)+5},2);
        bc_nodeset_temp = num_start:num_gap:num_end;
    else
        count2 = location_boundary(count)+2;
        while str2num(all{count2})  % Լ���ĵ㼯����������
            bc_nodeset_temp(flag) = any2num(all{count2},2);
            flag=flag+1;
            count2=count2+1;
        end
    end
    switch bctype{count}   % Լ�����Ͷ�Ӧ�����ɶ� ����ABAQUSԼ������������ʾ���õı߽���������
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
        case 'dfbdof'   %�������Լ������Լ�����͵Ĺؼ��֣����������Լ�������ɶȽ��ж�ȡ�������������������Լ����
            while str2num(all{location_boundary_temp(count)})   %����location_boundary_temp�������ԭ���ǿ��ܶ��Լ�����ǰ��յڶ��ַ�ʽ������
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
            bctypedof(count2)=[];    % ɾ�������ڵ�Լ��������ֻ��һ���ڵ�3�����ɶȣ���Ҫɾȥ��456���ɶȵ�Լ��
        else
            count2=count2+1;
        end
    end
    for count2 = 1:length(bc_nodeset_temp)
        BC=[BC (bc_nodeset_temp(count2)-1)*ndpn+bctypedof];  % ��Լ�������ɶ�
        count2=count2+1;
    end
end
fprintf('The total number constrained DOFs is :%d\n',length(BC));


%% ��ȡ������Ϣ
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

%% ��ȡ�غ�����---ֻ֧����ȡ�����غ�
 keyword_load={}; %�������¶���keyword_load -- �غɵĵ㼯����
for count = 1:length(location_load)
    keyword_load{count}=strcat('nset=',all{location_load(count)+1});
end
location_load_temp = location_load+1; % �غɵ㼯���ƴ���λ��  +1��λʩ���غɵ����ɶȵ�λ��
flag_load = 1;
for count = location_split:length(all)
    for count2=1:length(keyword_load)
        if strcmp(keyword_load{count2},all{count})
            location_load(flag_load) = count;  % �غɵ㼯��λ��
            flag_load = flag_load+1;
        end
    end
end

flag=1;
for count = 1:length(location_load)
    location_temp=location_load(count)+2;
    while any2num(all{location_temp},2)
        load_dof = any2num(all{location_load_temp(count)+1},1);%�غ�ʩ�ӵ����ɶ�
        load_node = any2num(all{location_temp},1);%�غɽڵ�
        load(flag,1)= (load_node-1)*ndpn+load_dof;
        load(flag,2) = any2num(all{location_load_temp(count)+2},2);
        location_temp=location_temp+1;
        flag=flag+1;
    end
end
end