%%%读取inp文件模型信息函数
function [Nodes, Elements] = Readmessage( fname )
fid = fopen(fname,'rt');  %fname文件名   r读取  t以txt格式打开
S = textscan(fid,'%s','Delimiter','\n'); %已经打开的文件  字符向量形式读取   'Delimiter','\n'分隔符为换行符  默认分隔符为空格
S = S{1};
%找到Node关键字所在的位置
idxS = strfind(S, 'Node');  %返回元胞数组 若数组中没有相应元素，则返回空
idx1 = find(not(cellfun(@isempty, idxS))); %cellfun(fun,A) 对元胞数组A分别使用函数fun   isempty(A) A为空返回逻辑值1   find 寻找非0元素的索引
%找到Element关键字所在的位置
idxS = strfind(S, 'Element');
idx2 = find(not(cellfun(@isempty, idxS)));
%找到Nset关键字所在位置
idxS = strfind(S, 'Nset');
idx3 = find(not(cellfun(@isempty, idxS)));
% 取出节点信息(元胞数组)
Nodes = S(idx1(1)+1:idx2(1)-1);  %以元胞数组形式取出
%将元胞数组转换为矩阵
Nodes = cell2mat(cellfun(@str2num,Nodes,'UniformOutput',false));  %'UniformOutput',false  以元胞形式返回输出值
% 取出单元(元胞数组)
elements = S(idx2+1:idx3(1)-1) ;
% 将元胞数组转换为矩阵
Elements = cell2mat(cellfun(@str2num,elements,'UniformOutput',false));
Nodes=Nodes(:,2:end);
Elements=Elements(:,2:end);
end