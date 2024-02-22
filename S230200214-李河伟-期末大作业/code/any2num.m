function [ number ] = any2num(any,index)
% index==1  --> int
% index==2  --> double
if index==1  % 提取any中的数字，忽略其他符号  12sy3,45,6 --> 123456
    number=str2num(any(isstrprop(any,'digit')));
end
if index==2
    if contains(any,'e')
        temp=strsplit(any,'e');
        a=str2double(regexp(temp{1}, '-?\d*\.?\d*', 'match'));
        b=str2double(regexp(temp{2}, '-?\d*\.?\d*', 'match'));
        number = a.*10.^b;
    else
        number=str2double(regexp(any, '-?\d*\.?\d*', 'match'));
    end
end
end