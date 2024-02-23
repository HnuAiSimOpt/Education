function [ff]=feasmbl_f(ff,f,index)
%-------------------------------------------------------------------
% Purpose:
%     Assembly of element force vector into the system force vector
%
% Synopsis:
%     [ff]=feasmbl_f(ff,f,index)
%
% Variable descriptions:
%     f - element force vector
%     ff - system force vector
%     index - d.o.f. vector associated with an element
%-------------------------------------------------------------------

edof=length(index);
for i=1:edof
    ii=index(i);
    ff(ii)=ff(ii)+f(i);
end

return