function dispcp(count,nes)
cpo = floor((count-1)/nes*100);
if cpo ~= floor(count/nes*100)
    cpo = floor(count/nes*100);
    fprintf('Computation Progress:%d%%\n',cpo);
end
end