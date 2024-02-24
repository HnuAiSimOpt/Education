function [selection]=selectElement()
fprintf('计算单元类型：\n');
fprintf('CPS6--二阶三角形平面应变单元\n');
fprintf('CPS4--线性四边形形平面应变单元\n');
fprintf('C3D4--线性四面体单元\n');
fprintf('S4R --缩减积分壳单元\n');
selection=input('\n','s');
switch selection
    case 'CPS6'
        fprintf('选择单元:%s\n',selection);
    case 'CPS4'
        fprintf('选择单元:%s\n',selection);
    case 'C3D4'
        fprintf('选择单元:%s\n',selection);
    case 'S4R'
        fprintf('选择单元:%s\n',selection);
    otherwise
        fprintf('输入错误\n');
end
% LOAD DATA
selection=strcat(selection,'.inp');

return