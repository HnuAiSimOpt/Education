function [xyz_ver,ind_ele]=comsol2matlab(file,n)
f=fopen(file); 
data_comsol=textscan(f,'%s','Delimiter','\n');
%% 获取节点
index_m=find(strcmp('# Mesh vertex coordinates',data_comsol{1,1}));
%网格节点数
n_ver=str2double(cell2mat(regexp(data_comsol{1,1}{index_m-3,1},'\-?\d*\.?\d*','match')));
%节点坐标
xyz_ver=str2num(char(data_comsol{1,1}{index_m+1:index_m+n_ver,1}));
%% 获取四面体单元
if n==10
    index_e=find(strcmp('10 # number of vertices per element',data_comsol{1,1}));
else
    index_e=find(strcmp('4 # number of vertices per element',data_comsol{1,1}));
end
%网格单元数
n_ver=str2double(cell2mat(regexp(data_comsol{1,1}{index_e+1,1},'\-?\d*\.?\d*','match')));

%单元节点索引
ind_ele=str2num(char(data_comsol{1,1}{index_e+3:index_e+n_ver+2,1}))+1;
end