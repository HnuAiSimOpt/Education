function [selection]=selectElement()
fprintf('���㵥Ԫ���ͣ�\n');
fprintf('CPS6--����������ƽ��Ӧ�䵥Ԫ\n');
selection=input('\n','s');
switch selection
    case 'CPS6'
        fprintf('ѡ��Ԫ:%s\n',selection);
    otherwise
        fprintf('�������\n');
end
% LOAD DATA
selection=strcat(selection,'.inp');

return