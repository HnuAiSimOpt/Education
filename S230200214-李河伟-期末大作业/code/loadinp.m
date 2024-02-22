
function [gcoords,nodes,nnpe,ndpn,ElasticM,PoissonR,BC,load,dimension,element_type,SecType,Thickness,Density] = loadinp(dictionary)
format long
%% Extract the info from INP exported from abaqus.
%      gcoords -- 节点坐标信息 [节点编号,节点坐标1,节点坐标2,节点坐标3...;]
%     ElasticM -- 弹性模量
%     PoissonR -- 泊松比
%        nodes -- 单元连接信息 [单元编号,节点1,节点2,节点3...;]
%           BC -- 边界条件（自由度）   默认边界条件约束为0
%         load -- 载荷 [承受载荷的自由度,载荷大小;]
%         nnpe -- 单元节点数，与单元类型有关
%         ndpn -- 节点自由度数
% all是以空格为分隔的读取字符，将所有的信息读取为一个细胞数组放在all中，通过all{?}访问

%% 文件准备---读取
fid_inp= fopen(dictionary,'r');
temp_all = textscan( fid_inp ,'%s');
all=temp_all{1};    % 解嵌套
clear temp_all;     % 释放内存
fclose(fid_inp);    % 关闭文件
clear fid_inp;
format long

%% 关键字及内存分配
keyword_node      = '*Node';         % node开始的标记——参考inp文件规范
keyword_element   = '*Element,';     % 单元开始
keyword_material  = '*Material,';    % 材料的关键字
keyword_load      = '*Cload';        % 集中载荷关键字
keyword_boundary  = '*Boundary';     % 约束的关键字
keyword_sec       = 'Section:';
keyword_split1    = '*End';          % 将前面part与后面的约束之类的分开
keyword_split2    = 'Part';
location_split    = 0;
location_node     = 0;
location_sec      = 0;
location_element  = 0;
location_material = 0;
location_load     = [];
location_boundary = [];
bc_nodeset_temp   = [];              % 约束点集
bctypedof         = [];              % 约束类型对应约束的自由度
flag_load         = 0;               % 可能有多种集中载荷
flag_boundary     = 0;               % 可能有多个约束条件
number_node       = 0;               % inp文件中记录的节点数
number_element    = 0;               % inp文件中记录的单元数，用于读取到的单元数进行核对
BC                = [];              % 被约束的自由度
bcset_smrk        = {};              %  
bctype            = [];              % 约束类型，支持abaqus约定的全部类型
gcoords           = zeros(1);
nodes             = zeros(1);
load              = zeros(1);        % [载荷所在自由度1,载荷大小1;
                                     %  载荷所在自由度2,载荷大小2;]
Thickness         = 0;               % SecType==solid --> Thickness=0
SecType           = {};
count             = 1;
Density           = 0;

%% 查找关键字的位置
while location_split==0
    if strcmp(keyword_split1,all{count})&&strcmp(keyword_split2,all{count+1}) % 将文件分断
        location_split = count;
    end
    count=count+1;
end
for count = 1:location_split
    if strcmp(keyword_node,all{count})&&location_node==0  %查找node开始的位置
        location_node = count;   %节点坐标开始的位置
    end
    if strcmp(keyword_element,all{count})&&location_element==0
        location_element = count;     %单元连接信息开始的地方
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
        location_load(flag_load+1) = count;   % 多个集中载荷或在一个节点同时施加多个集中载荷
        while ~strcmp('**',all{count2+4})
            flag_load = flag_load+1;
            location_load(flag_load+1) = count2+3;
            count2=count2+3;
        end
        flag_load=flag_load+1;
    end
    if strcmp(keyword_boundary,all{count})
        location_boundary(flag_boundary+1) = count;   % 多个约束
        flag_boundary = flag_boundary+1;
    end
end

%% 单元类型---根据单元类型确定单元节点数和自由度，进一步才能读取
element_type = all{location_element+1};
temp=strsplit(all{location_element+1},'=');
element_type=temp{2};
fprintf('Element type:%s\n',element_type);
switch element_type
    case 'S4R'             % 4节点四边形缩减积分壳单元
        nnpe=4;            % 单元节点数
        ndpn=6;            % 节点自由度数
        dimension=3;       % 三维单元 -- 确定每个节点的坐标数，xy or xyz
    case 'S4R5'            % 4节点四边形缩减积分薄壳单元，沙漏控制
        nnpe=4;
        ndpn=5;
        dimension=3;
    case 'type=S3R'        % 3节点三角形缩减积分壳单元
        nnpe=3;
        ndpn=6;
        dimension=3;
    case 'S4'              % 4节点通用壳单元 有限膜应变
        nnpe=4;
        ndpn=6;
        dimension=3;
    case 'S3'              % 3节点三角形壳单元 有限膜应变
        nnpe=3;
        ndpn=6;
        dimension=3;
    case 'C3D4'            % 4节点四面体线性单元
        nnpe=4;
        ndpn=3;
        dimension=3;
    case 'C3D6'            % 6节点线性三角形棱柱单元
        nnpe=6;
        ndpn=3;
        dimension=3;
    case 'C3D8'            % 8节点线性实体单元
        nnpe=8;
        ndpn=3;
        dimension=3;
    case 'C3D8R'           % 8节点线性实体单元 缩减积分 沙漏控制
        nnpe=8;
        ndpn=3;
        dimension=3;
    case 'C3D10M'          % 10节点修正四面体二次单元 沙漏控制
        nnpe=10;
        ndpn=3;
        dimension=3;
    case 'CPS6'            % 6节点三角形二次平面应力单元
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

%% 提取节点坐标参数
count = location_node+1;    % 节点坐标的起始位置是关键字的下一位
while count<location_element  % 节点坐标之后紧接着是单元的关键字
    column = rem(count-location_node,dimension+1);  % 列，当前位置对维度数加一取余，节点编号-节点坐标1-节点坐标2-节点坐标3
    row = floor((count-location_node)/(dimension+1))+1; % 行 向下取整后+1
    if column == 0  % 一行的最后一位的余数为0
        row = row-1;
        column = dimension+1;
    end
    gcoords(row,column) = any2num(all{count},2);   %提取小数
    count = count+1;
end

%% 提取单元节点连接信息
count = location_element+2;% 单元关键字 - 单元类型 - 单元信息开始   关键字到单元信息开始要+2
while str2num(all{count})  % 单元编号以及连接的节点编号全为正整数，在单元信息结束后为字符串，str2num之后为false，节点坐标可能为0
    column = rem(count-location_element-1,nnpe+1);
    row = floor((count-location_element-1)/(nnpe+1))+1;
    if column==0
        row=row-1;
        column=nnpe+1;
    end
    nodes(row,column)=any2num(all{count},1);    %提取整数
    count = count+1;
end

%% 验证节点数目和单元数目提取正确
if strcmp('generate',all{count+2})  %这里的count指向单元结束后的第一个位置--单元数，即*Nset
    num_end = any2num(all{count+4},2); %这里的count与上面提取单元信息的count紧密相连
    num_start = any2num(all{count+3},2);
    num_gap = any2num(all{count+5},2);
    number_node = (num_end-num_start)/num_gap+1;
    count=count+6;   % 这里count+6即指向标注单元数的位置，即*Elset，count还与下面读取单元数相连接
else
    while str2num(all{count+2})  %这是单元不是顺序生成的情况，这里count一直在移动，最后停在*Elset
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

%% 查找材料参数——弹性材料——只有一种材料
if strcmp('*Density',all{location_material+2})  % 仅限于弹性材料
    Density = any2num(all{location_material+3},2);
    location_material=location_material+2;
end
if strcmp('*Elastic',all{location_material+2})  % 仅限于弹性材料
    ElasticM=any2num(all{location_material+3},2);   % 提取小数
    PoissonR=any2num(all{location_material+4},2);   % 提取小数
else
    fprintf('Err:Only one elastic material is supported!\n');
end


%% 提取边界条件   支持ABAQUS两种定义约束的方式，提取边界条件模块已完善
% 先找到边界条件开始的位置，找出定义的点集和约束类型，其中约束类型分为按照约定的关键字定义和直接给出自由度来定义
location_boundary_temp=[];
for count = 1:flag_boundary  %这里count指示约束的数目  这里是找到约束的点集名称bcset_smrk与其对应的约束类型bctype
    bcset_smrk{count} = strcat('nset=',all{location_boundary(count)+1});  % 边界约束点集
    if str2num(all{location_boundary(count)+2})   % 如果这个位置是数字的话，说明是第二种定义约束的方式
        bctype{count} = 'dfbdof';   % 第二种定义约束的类型，直接给出约束的自由度，这里是自己约定了一个关键字，相当于另一种约束类型
        location_boundary_temp(count) = location_boundary(count)+2; %记录下这里直接给出的约束的自由度的位置，下面location_boundary变成了约束的点集的位置了
    else
        bctype{count} = all{location_boundary(count)+2};% 约束类型  bctype必须与bcset_smrk,location_boundary等对应
    end
end

for count = 1:length(all)
    for count2 = 1:length(bcset_smrk)
        if strcmp(bcset_smrk{count2},all{count})    % 查找约束点集的位置
            location_boundary(count2)=count;       % 这里location_boundary变成了约束点集的位置
        end
        count2=count2+1;
    end
    count=count+1;
end
flag=1;
for count=1:length(bctype)   %length(bctype)==length(location_boundary)  查找约束的自由度
    if strcmp('generate',all{location_boundary(count)+2})   % 约束的点集（连续）
        num_start = any2num(all{location_boundary(count)+3},2);
        num_end = any2num(all{location_boundary(count)+4},2);
        num_gap = any2num(all{location_boundary(count)+5},2);
        bc_nodeset_temp = num_start:num_gap:num_end;
    else
        count2 = location_boundary(count)+2;
        while str2num(all{count2})  % 约束的点集（非连续）
            bc_nodeset_temp(flag) = any2num(all{count2},2);
            flag=flag+1;
            count2=count2+1;
        end
    end
    switch bctype{count}   % 约束类型对应的自由度 根据ABAQUS约定的名称来表示常用的边界条件类型
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
        case 'dfbdof'   %这里是自己定义的约束类型的关键字，所以下面对约束的自由度进行读取，不是像上面的是事先约定的
            while str2num(all{location_boundary_temp(count)})   %这里location_boundary_temp是数组的原因是可能多个约束都是按照第二种方式给出的
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
            bctypedof(count2)=[];    % 删除不存在的约束，比如只有一个节点3个自由度，则要删去对456自由度的约束
        else
            count2=count2+1;
        end
    end
    for count2 = 1:length(bc_nodeset_temp)
        BC=[BC (bc_nodeset_temp(count2)-1)*ndpn+bctypedof];  % 被约束的自由度
        count2=count2+1;
    end
end
fprintf('The total number constrained DOFs is :%d\n',length(BC));


%% 提取截面信息
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

%% 提取载荷条件---只支持提取集中载荷
 keyword_load={}; %这里重新定义keyword_load -- 载荷的点集名称
for count = 1:length(location_load)
    keyword_load{count}=strcat('nset=',all{location_load(count)+1});
end
location_load_temp = location_load+1; % 载荷点集名称处的位置  +1即位施加载荷的自由度的位置
flag_load = 1;
for count = location_split:length(all)
    for count2=1:length(keyword_load)
        if strcmp(keyword_load{count2},all{count})
            location_load(flag_load) = count;  % 载荷点集的位置
            flag_load = flag_load+1;
        end
    end
end

flag=1;
for count = 1:length(location_load)
    location_temp=location_load(count)+2;
    while any2num(all{location_temp},2)
        load_dof = any2num(all{location_load_temp(count)+1},1);%载荷施加的自由度
        load_node = any2num(all{location_temp},1);%载荷节点
        load(flag,1)= (load_node-1)*ndpn+load_dof;
        load(flag,2) = any2num(all{location_load_temp(count)+2},2);
        location_temp=location_temp+1;
        flag=flag+1;
    end
end
end