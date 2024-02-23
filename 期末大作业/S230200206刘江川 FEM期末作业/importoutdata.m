function comsolData = importoutdata(filename, dataLines)
%IMPORTFILE 从文本文件中导入数据
%  comsolData = IMPORTFILE(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。  返回数值数据。
%
%  comsolData = IMPORTFILE(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME
%  中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  comsolData = importfile("E:\Demo\Matlab\2023_FEM\位移和力.txt", [10, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2024-02-16 11:09:26 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
    dataLines = [10, Inf];
end

%% 设置导入选项并导入数据
opts = delimitedTextImportOptions("NumVariables", 10);

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = " ";

% 指定列名称和类型
opts.VariableNames = ["VarName1", "Model", "FEMmph", "VarName4", "VarName5", "VarName6", "VarName7", "Var8", "Var9", "Var10"];
opts.SelectedVariableNames = ["VarName1", "Model", "FEMmph", "VarName4", "VarName5", "VarName6", "VarName7"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "string", "string", "string"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% 指定变量属性
opts = setvaropts(opts, ["Var8", "Var9", "Var10"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var8", "Var9", "Var10"], "EmptyFieldRule", "auto");

% 导入数据
comsolData = readtable(filename, opts);

%% 转换为输出类型
comsolData = table2array(comsolData);
end