function [selection]=selectElement()
fprintf('计算单元类型：\n');
fprintf('CPS6--二阶三角形平面应变单元\n');
selection=input('\n','s');
switch selection
    case 'CPS6'
        fprintf('选择单元:%s\n',selection);
    otherwise
        fprintf('输入错误\n');
end
% LOAD DATA
selection=strcat(selection,'.inp');

return