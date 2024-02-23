%-----------------------------*单元长度*----------------------------------
function L = eLength(DYJD,JDZB)
ZDY = length(DYJD);
L = zeros(1,ZDY);
for i = 1:ZDY
JDZB1 = JDZB(DYJD(1,i),:);
JDZB2 = JDZB(DYJD(2,i),:);
L(i) = sqrt(sum(power(JDZB1 - JDZB2,2)));
end