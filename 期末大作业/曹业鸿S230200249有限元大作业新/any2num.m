function [ number ] = any2num(any,index)
if index==1
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